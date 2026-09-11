"""Regressions observed during native multi-module pilot validation."""
from io import StringIO
from argparse import Namespace

import pytest
from rich.cells import cell_len
from rich.console import Console

from postgwas.core.checkpointing import _package_content_identity
from postgwas.core import checkpointing
from postgwas.core.ui.progress import StageProgress


@pytest.mark.parametrize('suffix', ('.py', '.R', '.yaml', '.json'))
def test_package_identity_changes_with_source_and_bundled_resources(tmp_path, suffix):
    resource = tmp_path / ('analysis' + suffix)
    resource.write_text('first\n')
    original = _package_content_identity(tmp_path)
    resource.write_text('other\n')  # same file size; content, not mtime, determines identity
    assert _package_content_identity(tmp_path) != original
    resource.write_text('first\n')
    assert _package_content_identity(tmp_path) == original
    resource.rename(tmp_path / ('renamed' + suffix))
    assert _package_content_identity(tmp_path) != original


def test_package_identity_ignores_python_cache_but_tracks_empty_files(tmp_path):
    (tmp_path / '__init__.py').touch()
    original = _package_content_identity(tmp_path)
    cache = tmp_path / '__pycache__'
    cache.mkdir()
    (cache / 'module.cpython-310.pyc').write_bytes(b'bytecode')
    assert _package_content_identity(tmp_path) == original
    (tmp_path / 'model.json').write_text('{}')
    assert _package_content_identity(tmp_path) != original


def test_package_identity_rejects_empty_package(tmp_path):
    with pytest.raises(RuntimeError, match='Cannot fingerprint'):
        _package_content_identity(tmp_path)


def test_software_identity_includes_live_package_content(monkeypatch):
    content = {"files": 1, "sha256": "first"}
    monkeypatch.setattr(checkpointing, '_package_content_identity', lambda root: dict(content))
    original = checkpointing.software_identity()
    assert original['postgwas_content'] == content
    content['sha256'] = 'changed'
    assert checkpointing.software_identity() != original


@pytest.mark.parametrize('layout_key', ('output_layout', 'fine_mapping_output_layout'))
@pytest.mark.parametrize('namespace', (False, True))
def test_input_discovery_never_expands_output_layouts(tmp_path, monkeypatch, layout_key, namespace):
    monkeypatch.chdir(tmp_path)
    output = tmp_path / 'results'
    output.mkdir()
    live_log = output / 'screen.log'
    live_log.write_text('active\n')
    upstream = output / 'upstream.tsv'
    upstream.write_text('validated upstream data\n')
    fields = {layout_key: {'results_directory': 'results'}, 'input_file': upstream}
    value = Namespace(**fields) if namespace else {'resolved_configuration': fields}
    assert set(checkpointing.discover_input_files(value).values()) == {upstream}
    live_log.write_text('another run appended\n')
    assert set(checkpointing.discover_input_files(value).values()) == {upstream}


def test_input_discovery_prunes_output_directory_but_keeps_reference_directory(tmp_path):
    output = tmp_path / 'outputs'
    output.mkdir()
    (output / 'live.log').write_text('active\n')
    reference = tmp_path / 'reference'
    reference.mkdir()
    resource = reference / 'data.tsv'
    resource.write_text('resource\n')
    value = Namespace(fine_mapping_output_directory=output, reference_directory=reference)
    assert set(checkpointing.discover_input_files(value).values()) == {resource}


@pytest.mark.parametrize('width', (80, 120))
@pytest.mark.parametrize('nested', (False, True))
def test_stage_outcome_long_values_keep_continuation_alignment(width, nested):
    stream = StringIO()
    progress = StageProgress('Test', enabled=True, outcome_label_width=42,
                             console=Console(file=stream, width=width, color_system=None))
    value = 'ADHD2022_iPSYCH_deCODE_PGC_long_input_filename.vcf.gz'
    fields = [('analysis', 'Input')] if nested else []
    fields.append(('info', 'Input VCF', value))
    progress.start_step(1, 1, 'Read input')
    progress.complete_step(1, 1, 'Read input', outcome_fields=fields)
    lines = stream.getvalue().splitlines()
    first = next(i for i, line in enumerate(lines) if 'Input VCF' in line)
    prefix, fragment = lines[first].split(' : ', 1)
    value_column = cell_len(prefix + ' : ')
    pieces = [fragment]
    for line in lines[first + 1:]:
        if 'All 1 stages completed' in line:
            break
        assert line.startswith(' ' * value_column), repr(line)
        pieces.append(line[value_column:])
        assert cell_len(line) <= width
    assert ''.join(pieces) == value
    assert cell_len(lines[first]) <= width
