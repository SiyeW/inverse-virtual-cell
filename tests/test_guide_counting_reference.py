import csv
import hashlib
import io
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import prepare_guide_counting_reference as prep


def fixture(tmp_path):
    root = tmp_path/'source';root.mkdir()
    path = root/'features.csv'
    rows = [dict(zip(prep.CAPTURE+prep.TARGET, values)) for values in [
        ['PB.x-1','PB.x-1','R2','ACGT(BC)','ACGT','CRISPR Guide Capture','PB.x','PB.x'],
        ['NTC-1','NTC-1','R2','ACGT(BC)','TGCA','CRISPR Guide Capture','Non-Targeting','Non-Targeting']]]
    path.write_bytes(prep.encode(rows, prep.CAPTURE+prep.TARGET))
    return path, rows


def run(tmp_path, source):
    return prep.prepare(source, hashlib.sha256(source.read_bytes()).hexdigest(), tmp_path/'out', tmp_path)


def test_capture_and_authority_remain_exact(tmp_path):
    source, rows = fixture(tmp_path);before = source.read_bytes()
    receipt = run(tmp_path, source)
    out = tmp_path/'out'
    assert receipt['guides'] == 2 and receipt['target_and_control_relabeling'] is False
    assert source.read_bytes() == before == (out/'source_feature_reference.csv').read_bytes()
    capture = list(csv.DictReader(io.StringIO((out/'guide_counting_features.csv').read_text())))
    assert capture == [{k:r[k] for k in prep.CAPTURE} for r in rows]
    mapping = list(csv.DictReader(io.StringIO((out/'authoritative_target_map.csv').read_text())))
    assert mapping[0]['target_gene_id'] == 'PB.x'
    assert mapping[1]['target_gene_id'] == 'Non-Targeting'
    for name, proof in receipt['files'].items():
        assert hashlib.sha256((out/name).read_bytes()).hexdigest() == proof['sha256']
    with pytest.raises(FileExistsError):run(tmp_path, source)


def test_bad_source_hash_and_source_overlap_stop(tmp_path):
    source, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match='SHA256 differs'):
        prep.prepare(source, '0'*64, tmp_path/'out', tmp_path)
    with pytest.raises(ValueError, match='overlaps'):
        prep.prepare(source, hashlib.sha256(source.read_bytes()).hexdigest(), source.parent/'out', tmp_path)
    assert not (tmp_path/'out').exists()


def test_relative_paths_are_rejected(tmp_path):
    source, _ = fixture(tmp_path)
    with pytest.raises(ValueError, match='Absolute'):
        prep.prepare(source, hashlib.sha256(source.read_bytes()).hexdigest(), 'relative-output', tmp_path)


@pytest.mark.parametrize('change', ['duplicate', 'extra_column', 'bad_sequence', 'wrong_type', 'missing_value'])
def test_invalid_inputs_fail_before_output(tmp_path, change):
    source, rows = fixture(tmp_path)
    fields = prep.CAPTURE+prep.TARGET
    if change == 'duplicate':rows[1]['id'] = rows[0]['id']
    if change == 'extra_column':
        fields = fields + ['unknown']
        for r in rows:r['unknown'] = 'do not silently discard'
    if change == 'bad_sequence':rows[0]['sequence'] = '???'
    if change == 'wrong_type':rows[0]['feature_type'] = 'Antibody Capture'
    if change == 'missing_value':rows[0]['pattern'] = ''
    source.write_bytes(prep.encode(rows, fields))
    with pytest.raises(ValueError):run(tmp_path, source)
    assert not (tmp_path/'out').exists()
