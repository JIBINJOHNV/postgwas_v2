"""Whole-study replacement evidence uses raw study AF and original-EA replacement AF."""
import math

import polars as pl
import pytest

from postgwas.modules.harmonisation import allele_frequency as af
from postgwas.modules.harmonisation.policies import default_policies


def compare(a, b, *, total=None, overrides=None, batches=None, **kwargs):
    frame = pl.DataFrame({'A': a, 'B': b}, schema={'A': pl.Float64, 'B': pl.Float64})
    paired = frame if batches is None else lambda: frame.iter_slices(batches)
    return af.validate_replacement_frequency_pairs(paired, 'A', 'B',
        default_policies().with_overrides(overrides or {}),
        study_variant_count=len(a) if total is None else total, **kwargs)


def variable(n=1000):
    return [.1 + .3 * i / max(n - 1, 1) for i in range(n)]


def test_replacement_folds_only_b_and_preserves_both_columns():
    frame = pl.DataFrame({'A': [.2, .3, .4], 'B': [.79, .69, .59]})
    before = frame.clone()
    evidence = af.validate_replacement_frequency_pairs(frame, 'A', 'B',
        default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2}),
        study_variant_count=3, source={'file': '/reference.tsv', 'column': 'ALT_AF', 'mapping': {'a1': 'ALT'}})
    assert evidence['status'] == 'passed'
    assert evidence['mean_absolute_difference'] == pytest.approx(.01)
    assert evidence['correlation'] == pytest.approx(1.)
    assert evidence['major_count'] == 3 and evidence['major_fraction'] == 1.
    assert evidence['source']['mapping'] == {'a1': 'ALT'}
    assert frame.equals(before)
    raw = compare([.8, .3, .4], [.8, .7, .6],
        overrides={'eaf.maf_reference_min_overlap': 2})
    assert raw['status'] == 'inconclusive'
    assert raw['correlation'] < 0


@pytest.mark.parametrize('count,total,want,reason', [
    (999, 999, 'inconclusive', 'insufficient_overlap'),
    (1000, 2000, 'passed', None),
    (1000, 2001, 'inconclusive', 'insufficient_overlap_fraction'),
    (0, 1, 'inconclusive', 'insufficient_overlap'),
])
def test_replacement_count_and_fixed_coverage_boundaries(count, total, want, reason):
    a = variable(count)
    evidence = compare(a, [1 - v for v in a], total=total)
    assert evidence['status'] == want
    assert evidence['comparable_variants'] == count
    assert evidence['study_variant_count'] == total
    assert evidence['observed_overlap_fraction'] == count / total
    assert evidence['minimum_overlap'] == 1000
    assert evidence['minimum_overlap_fraction'] == .5
    if reason:
        assert reason in evidence['failed_requirements']
        assert evidence['correlation'] is None
        assert evidence['mean_absolute_difference'] is None


@pytest.mark.parametrize('major,want', [(199, 'inconclusive'), (200, 'inconclusive'), (201, 'passed')])
def test_replacement_major_fraction_is_strict_and_half_is_not_major(major, want):
    a = variable(999) + [.5]
    b = [1 - v if i < major else v for i, v in enumerate(a)]
    evidence = compare(a, b)
    assert evidence['major_count'] == major
    assert evidence['major_fraction'] == major / 1000
    assert evidence['minimum_major_fraction'] == .20
    assert evidence['status'] == want
    assert evidence['failed_requirements'] == ([] if want == 'passed' else ['major_fraction_not_above_minimum'])


@pytest.mark.parametrize('delta,want', [(0., 'passed'), (1e-9, 'inconclusive')])
def test_replacement_correlation_threshold_equality(delta, want):
    a = variable()
    b = [1 - (v + .02 * math.sin(i)) for i, v in enumerate(a)]
    expected = pl.DataFrame({'a': a, 'c': [1 - v for v in b]}).select(pl.corr('a', 'c')).item()
    threshold = compare(a, b)['correlation']
    evidence = compare(a, b, overrides={'eaf.study_reference_min_correlation': threshold + delta})
    assert evidence['correlation'] == pytest.approx(expected, abs=1e-14)
    assert evidence['status'] == want


@pytest.mark.parametrize('threshold,want', [(.125, 'passed'), (.125 - 1e-9, 'inconclusive')])
def test_replacement_mean_absolute_difference_threshold_equality(threshold, want):
    evidence = compare([.25, .5] * 500, [.875, .625] * 500,
        overrides={'eaf.maf_reference_max_mean_absolute_difference': threshold})
    assert evidence['mean_absolute_difference'] == .125
    assert evidence['status'] == want


@pytest.mark.parametrize('method', ['pearson', 'spearman'])
@pytest.mark.parametrize('a,b', [([.2] * 1000, [.8] * 1000),
    (variable(), [.25, .75] * 500)])
def test_replacement_constant_vectors_do_not_supply_correlation(method, a, b):
    evidence = compare(a, b, overrides={'eaf.maf_reference_correlation_method': method})
    assert evidence['correlation'] is None
    assert evidence['status'] == 'inconclusive'
    assert evidence['failed_requirements'][0] == 'correlation_unavailable'


def test_replacement_missing_invalid_and_diagnostic_exclusions_share_n():
    frame = pl.DataFrame({'A': [.1, .2, .3, None, .2, float('nan'), .2, -1., .2],
        'B': [.9, .8, .7, .8, None, .8, float('inf'), .8, 1.1],
        'eligible': [True] * 8 + [False]})
    evidence = af.validate_replacement_frequency_pairs(frame, 'A', 'B',
        default_policies().with_overrides({'eaf.maf_reference_min_overlap': 2,
            'eaf.maf_reference_min_overlap_fraction': .3}), study_variant_count=10,
        eligibility_column='eligible', excluded_counts={'palindromic_rows': 1, 'non_biallelic_rows': 1})
    assert evidence['comparable_variants'] == 3
    assert evidence['study_variant_count'] == 10
    assert evidence['observed_overlap_fraction'] == .3
    assert evidence['excluded_counts'] == {'palindromic_rows': 1, 'non_biallelic_rows': 1,
        'diagnostic_ineligible_rows': 1, 'missing_frequency_rows': 2,
        'non_finite_frequency_rows': 2, 'out_of_range_frequency_rows': 1}
    assert evidence['major_count'] == 3 and evidence['major_fraction'] == 1.
    assert evidence['status'] == 'passed'


def test_replacement_reports_all_measured_failed_requirements_without_classifier(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('original dual/legacy classifier must not run')
    monkeypatch.setattr(af, '_compare_frequency_pairs', forbidden)
    evidence = compare([.99] * 1000, [.1] * 1000,
        overrides={'eaf.maf_reference_error_margin': .5})
    assert evidence['failed_requirements'] == ['correlation_unavailable',
        'mean_absolute_difference_above_maximum', 'major_fraction_not_above_minimum']
    assert evidence['status'] == 'inconclusive'


@pytest.mark.parametrize('total', [0, -1, True, None, 1.5])
def test_replacement_rejects_invalid_fixed_denominator(total):
    frame = pl.DataFrame({'A': [.1, .2], 'B': [.9, .8]})
    with pytest.raises(af.StudyAlleleFrequencyError, match='positive integer'):
        af.validate_replacement_frequency_pairs(frame, 'A', 'B', study_variant_count=total)


def test_replacement_rejects_more_input_pairs_than_study_rows():
    with pytest.raises(af.StudyAlleleFrequencyError, match='exceed'):
        compare([.1, .2, .3], [.9, .8, .7], total=2)


@pytest.mark.parametrize('batch_size', [1, 17, 333, 1000])
def test_chunked_pearson_pools_centered_statistics_not_correlations(batch_size):
    a = [.2 + 1e-7 * i for i in range(1000)]
    b = [1 - (v + 1e-8 * math.sin(i)) for i, v in enumerate(a)]
    full = compare(a, b)
    chunked = compare(a, b, batches=batch_size)
    assert chunked['correlation'] == pytest.approx(full['correlation'], abs=1e-11)
    assert chunked['mean_absolute_difference'] == pytest.approx(full['mean_absolute_difference'], abs=1e-15)
    assert chunked['major_fraction'] == full['major_fraction']
    assert chunked['status'] == full['status'] == 'passed'


def test_chunked_spearman_uses_one_global_rank_with_ties():
    a = [.1, .4, .2, .2, .3, .1] * 200
    b = [.6, .9, .7, .7, .8, .6] * 200
    expected = pl.DataFrame({'a': a, 'c': [1 - v for v in b]}).select(pl.corr('a', 'c', method='spearman')).item()
    evidence = compare(a, b, batches=7,
        overrides={'eaf.maf_reference_correlation_method': 'spearman'}, spearman_memory_budget_bytes=16000000)
    assert evidence['correlation_method'] == 'spearman'
    assert evidence['correlation'] == pytest.approx(expected)
    assert evidence['correlation'] == pytest.approx(-1.)


def test_chunked_spearman_rejects_insufficient_resource_budget():
    with pytest.raises(af.StudyAlleleFrequencyError, match='Spearman.*memory budget'):
        compare(variable(), [.8] * 1000, batches=10, spearman_memory_budget_bytes=100,
            overrides={'eaf.maf_reference_correlation_method': 'spearman'})


def test_count_failure_does_not_rank_or_measure(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError('no metric or rank work before evidence count gate')
    monkeypatch.setattr(af, '_measure_frequency_pairs', forbidden, raising=False)
    evidence = compare([.2], [.8], total=100)
    assert evidence['failed_requirements'] == ['insufficient_overlap', 'insufficient_overlap_fraction']


def test_spearman_budget_is_checked_before_replay():
    calls = []
    frame = pl.DataFrame({'A': variable(), 'B': [1 - v for v in variable()]})
    def batches():
        calls.append('read')
        assert len(calls) == 1, 'budget check must precede second-pass materialization'
        yield frame
    with pytest.raises(af.StudyAlleleFrequencyError, match='Spearman.*memory budget'):
        af.validate_replacement_frequency_pairs(batches, 'A', 'B',
            default_policies().with_overrides({'eaf.maf_reference_correlation_method': 'spearman',
                'execution.memory_budget_gb': .000001}), study_variant_count=1000)


def test_spearman_direct_library_null_budget_preserves_exact_global_method():
    evidence = compare(variable(), [1 - v for v in variable()], batches=17,
        overrides={'eaf.maf_reference_correlation_method': 'spearman'})
    assert evidence['correlation'] == pytest.approx(1.)
    assert evidence['status'] == 'passed'


def test_chunked_metrics_are_invariant_to_partition_order():
    a = variable()
    b = [1 - (v + .01 * math.sin(i)) for i, v in enumerate(a)]
    frame = pl.DataFrame({'A': a, 'B': b})
    evidence = af.validate_replacement_frequency_pairs(lambda: reversed(list(frame.iter_slices(31))),
        'A', 'B', study_variant_count=1000)
    full = compare(a, b)
    for key in ('comparable_variants', 'major_count', 'major_fraction', 'excluded_counts', 'status'):
        assert evidence[key] == full[key]
    assert evidence['correlation'] == pytest.approx(full['correlation'], abs=1e-14)
    assert evidence['mean_absolute_difference'] == pytest.approx(full['mean_absolute_difference'], abs=1e-15)


def test_chunked_pearson_matches_high_precision_on_near_constant_data():
    from decimal import Decimal, localcontext
    a = [.499999 + 1e-10 * i for i in range(1000)]
    b = [1 - (v + 1e-11 * math.sin(i)) for i, v in enumerate(a)]
    with localcontext() as precision:
        precision.prec = 60
        x = [Decimal(v) for v in a]
        y = [Decimal(1 - v) for v in b]
        mx, my = sum(x) / 1000, sum(y) / 1000
        covariance = sum((v - mx) * (w - my) for v, w in zip(x, y))
        expected = float(covariance / (sum((v - mx)**2 for v in x) * sum((w - my)**2 for w in y)).sqrt())
    evidence = compare(a, b, batches=17)
    assert evidence['correlation'] == pytest.approx(expected, abs=1e-11)
    assert evidence['status'] == 'passed'


def test_tiny_nonconstant_valid_frequencies_have_finite_pearson():
    a = [(i + 1) * 1e-100 for i in range(1000)]
    evidence = compare(a, a, batches=17)
    assert evidence['correlation'] == pytest.approx(1.)
    assert evidence['failed_requirements'] == ['major_fraction_not_above_minimum']


@pytest.mark.parametrize('method', ['pearson', 'spearman'])
def test_changed_replay_count_is_invalid_comparison(method):
    frame = pl.DataFrame({'A': variable(), 'B': [1 - v for v in variable()]})
    calls = []
    def batches():
        calls.append('read')
        yield frame if len(calls) == 1 else frame.head(999)
    with pytest.raises(af.StudyAlleleFrequencyError, match='changed.*pair count'):
        af.validate_replacement_frequency_pairs(batches, 'A', 'B',
            default_policies().with_overrides({'eaf.maf_reference_correlation_method': method}),
            study_variant_count=1000)


def test_replacement_evidence_and_thresholds_reach_canonical_log(tmp_path):
    from postgwas.core.pipeline_logging import PipelineLogger
    with PipelineLogger('harmonisation', 'replacement', str(tmp_path / 'logs')) as logger:
        evidence = compare(variable(), [1 - v for v in variable()], logger=logger,
            source={'file': '/declared/reference.tsv', 'column': 'ALT_AF'})
    text = next((tmp_path / 'logs').glob('*replacement.log')).read_text()
    for key in ('Whole-study replacement AF validation', 'minimum_overlap', 'minimum_overlap_fraction',
            'minimum_correlation', 'maximum_mean_absolute_difference', 'minimum_major_fraction',
            'major_count', 'failed_requirements', '/declared/reference.tsv', 'ALT_AF'):
        assert key in text
    assert evidence['status'] == 'passed'


def test_folded_agreement_does_not_prove_declared_allele_direction():
    a = [.1, .2, .3, .4] * 250
    b = [.1, .8, .3, .6] * 250
    for candidate in (b, [1 - value for value in b]):
        evidence = compare(a, candidate)
        assert evidence['status'] == 'passed'
        assert evidence['correlation'] == pytest.approx(1., abs=1e-14)
        assert evidence['mean_absolute_difference'] == pytest.approx(0., abs=1e-15)
        assert evidence['major_fraction'] == .5


def test_unequal_pools_do_not_average_perfect_local_correlations():
    a = [.05, .10] * 5 + [.30, .35] * 495
    c = [.30, .35] * 5 + [.05, .10] * 495
    frame = pl.DataFrame({'A': a, 'B': [1 - value for value in c]})
    policies = default_policies().with_overrides({'eaf.maf_reference_max_mean_absolute_difference': .5})
    for paired in (frame, lambda: iter((frame.head(10), frame.tail(990)))):
        evidence = af.validate_replacement_frequency_pairs(paired, 'A', 'B',
            policies, study_variant_count=1000)
        # Hand-derived covariance .00000625 / variance .00124375 = 1/199;
        # the two separate chromosome correlations would each equal one.
        assert evidence['correlation'] == pytest.approx(1 / 199, abs=1e-14)
        assert evidence['mean_absolute_difference'] == pytest.approx(.25, abs=1e-15)
        assert evidence['failed_requirements'] == ['correlation_below_minimum']
        assert evidence['status'] == 'inconclusive'
