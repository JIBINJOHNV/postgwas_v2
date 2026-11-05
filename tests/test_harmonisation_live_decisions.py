"""Completed study decisions are visible before later preparation work starts."""

import pytest

from postgwas.modules.harmonisation import service
from test_harmonisation_study_report_integration import _prepare


def test_settled_study_decisions_are_visible_inside_partition_writer(tmp_path, monkeypatch, capsys):
    write = service.write_chromosome_partitions
    observed = []
    def checked_write(*args, **kwargs):
        text = capsys.readouterr().out
        observed.append(text)
        assert 'Selected build' in text
        assert 'Selected consensus' in text
        assert 'Selected type' in text
        assert 'Declared study EAF retained' in text
        return write(*args, **kwargs)
    monkeypatch.setattr(service, 'write_chromosome_partitions', checked_write)
    prepared = _prepare(tmp_path, monkeypatch)
    tail = capsys.readouterr().out
    full = ''.join(observed) + tail
    assert full.count('Study-wide decisions and evidence') == 1
    assert full.count('Selected build') == 1
    assert full.count('Selected consensus') == 1
    assert full.count('🧬  Allele frequency') == 1
    assert full.count('Imputation-quality score') == 1
    assert prepared.study_decisions['eaf_reference_validation']['status'] == 'not_required'


def test_late_failure_does_not_replay_published_decisions(tmp_path, monkeypatch, capsys):
    with pytest.raises(OSError, match='partition fixture failure'):
        _prepare(tmp_path, monkeypatch, failure_stage='partition')
    text = capsys.readouterr().out
    assert text.count('Selected build') == 1
    log = (tmp_path / 'logs' / 'study_dataset.log').read_text()
    assert log.index('Selected build') < log.index('Prepare chromosome files')
    assert text.count('Study-wide decisions and evidence') == 1
    assert 'Dataset stopped' in text
    assert 'partition fixture failure' in text
