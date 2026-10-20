import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[1]/'scripts'/'audit_finalized_delivery.py'
spec = importlib.util.spec_from_file_location('independent_auditor', SCRIPT)
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def write_json(path, value):
    path.write_text(json.dumps(value))
    return path.read_bytes()


def fixture(tmp_path, monkeypatch):
    source, target, control = [tmp_path/name for name in ('source', 'target', 'control')]
    for p in (source, target, control):
        p.mkdir()
    for p in (source/'cellrangerresult', target/'cellrangerresult'):
        p.mkdir()
        (p/'tiny').write_bytes(b'abc')
    src = a.regular_identity(source/'cellrangerresult/tiny')
    dst = a.regular_identity(target/'cellrangerresult/tiny')
    vendor = hashlib.md5(b'abc').hexdigest()
    sha = hashlib.sha256(b'abc').hexdigest()
    rel = 'cellrangerresult/tiny'
    manifest = dict(source=str(source), target=str(target), control=str(control),
                    cleanup_sources=False, files={rel: src}, directories=['cellrangerresult'],
                    total_files=1, total_bytes=3, vendor_md5={rel: vendor})
    manifest_sha = a.digest(write_json(control/'manifest.json', manifest))
    original_identity = dict(dst, inode=dst['inode']+1)
    old = dict(source=dict(identity=src, sha256=sha, md5=vendor),
               target=dict(identity=original_identity, sha256=sha, md5=vendor),
               vendor_md5=vendor, accepted_at='2026-10-20T00:00:00+08:00')
    prior = dict(manifest_sha256=manifest_sha, verified={rel: old})
    prior_raw = json.dumps(prior).encode()
    prior_sha = a.digest(prior_raw)
    operation = control/'handoff'
    before = operation/'before'
    before.mkdir(parents=True)
    (before/'ledger.json').write_bytes(prior_raw)
    owner = dict(pid=123, started_at='before', manifest_sha256=manifest_sha)
    owner_raw = write_json(before/'owner.before.json', owner)
    error_raw = b'Traceback (most recent call last):\nKeyboardInterrupt\n'
    (before/'run.stderr.log').write_bytes(error_raw)
    snapshots = {n: dict(bytes=len(raw), sha256=a.digest(raw)) for n, raw in (
        ('ledger.json', prior_raw), ('owner.before.json', owner_raw), ('run.stderr.log', error_raw))}
    write_json(before/'snapshot.json', snapshots)
    gate = dict(stage='old_writer_exited', manifest_sha256=manifest_sha,
                content_receipt_files=1, ledger_before_signal_sha256=prior_sha,
                previous_pid=123, previous_started_at='before', signalled=True)
    write_json(operation/'status.json', gate)
    history = control/'finalize_history'/(prior_sha+'_1')
    history.mkdir(parents=True)
    (history/'ledger.json').write_bytes(prior_raw)
    new = copy.deepcopy(old)
    new['target']['identity'] = dst
    new['target_content_identity'] = original_identity
    new['target_inode_audits'] = [dict(previous=original_identity, current=dst,
        time='now', policy=a.POLICY, sha256=sha, md5=vendor,
        method='inode-only metadata audit; existing full content hashes retained')]
    ledger = dict(manifest_sha256=manifest_sha, target_identity_policy=a.POLICY, verified={rel: new})
    ledger_sha = a.digest(write_json(control/'ledger.json', ledger))
    acceptance = dict(complete=True, finalize_only=True, source_retained=True,
        data_files_rehashed=0, total_files=1, total_bytes=3, source=str(source), target=str(target),
        vendor_md5_passed=1, manifest_sha256=manifest_sha, ledger_sha256=ledger_sha,
        target_identity_policy=a.POLICY, prior_ledger_sha256=prior_sha,
        inode_audited_files=1, inode_audited_bytes=3)
    write_json(control/'final_acceptance.json', acceptance)
    write_json(control/'status.json', dict(stage='complete', finalize_only=True))
    monkeypatch.setattr(a, 'mount_type', lambda _: 'cifs')
    monkeypatch.setattr(a, 'process_identity', lambda _: None)
    original_read = Path.read_bytes
    def bounded_control_only(path):
        if source in path.parents or target in path.parents:
            pytest.fail('Independent audit must never read a data body')
        return original_read(path)
    monkeypatch.setattr(Path, 'read_bytes', bounded_control_only)
    return control, manifest_sha, source, target


def test_independent_audit_preserves_all_records_and_never_reads_data(tmp_path, monkeypatch):
    control, sha, _, _ = fixture(tmp_path, monkeypatch)
    before = {str(p): p.read_bytes() for p in control.rglob('*') if p.is_file()}
    result = a.audit(control, sha, 'handoff', [(123, 42)])
    assert result['complete'] and result['source_identity_exact']
    assert result['data_body_bytes_read'] == result['data_files_rehashed'] == 0
    assert result['accepted_inode_audited_files'] == 1
    assert before == {str(p): p.read_bytes() for p in control.rglob('*') if p.is_file()}


@pytest.mark.parametrize('change', ['lock', 'partial', 'directory', 'source', 'target',
                                  'archive', 'ledger', 'live_actor', 'mount'])
def test_auditor_fails_closed_without_repairing_anything(tmp_path, monkeypatch, change):
    control, sha, source, target = fixture(tmp_path, monkeypatch)
    if change == 'lock': (control/'owner.lock').write_text('owned')
    elif change == 'partial': (target/'new.receiving').write_bytes(b'keep')
    elif change == 'directory': (target/'extra').mkdir()
    elif change == 'source': (source/'cellrangerresult/tiny').write_bytes(b'changed')
    elif change == 'target': (target/'cellrangerresult/tiny').write_bytes(b'changed')
    elif change == 'archive': (control/'handoff/before/ledger.json').write_text('{}')
    elif change == 'ledger': (control/'ledger.json').write_text('{}')
    elif change == 'live_actor':
        monkeypatch.setattr(a, 'process_identity', lambda _: dict(state='S', startticks=42))
    elif change == 'mount': monkeypatch.setattr(a, 'mount_type', lambda _: 'ext4')
    before = {str(p): p.read_bytes() for p in control.rglob('*') if p.is_file()}
    with pytest.raises((RuntimeError, KeyError)):
        a.audit(control, sha, 'handoff', [(123, 42)])
    assert before == {str(p): p.read_bytes() for p in control.rglob('*') if p.is_file()}


def test_post_acceptance_inode_remap_is_separate_observation(tmp_path, monkeypatch):
    control, sha, _, target = fixture(tmp_path, monkeypatch)
    original = a.regular_identity
    def remapped(path):
        value = original(path)
        if target in Path(path).parents:
            value['inode'] += 10
        return value
    monkeypatch.setattr(a, 'regular_identity', remapped)
    result = a.audit(control, sha, 'handoff')
    assert result['complete'] and len(result['post_acceptance_inode_only_observations']) == 1
    assert result['data_files_rehashed'] == 0


@pytest.mark.parametrize('change', ['content_hash', 'accepted_at', 'stable_field', 'chain'])
def test_audit_chain_rejects_modified_original_proof(change):
    identity = dict(bytes=3, mtime_ns=42, inode=1, device=2)
    old = dict(target=dict(identity=identity, sha256='a'*64, md5='b'*32), accepted_at='first')
    new = copy.deepcopy(old)
    if change == 'content_hash': new['target']['sha256'] = 'c'*64
    elif change == 'accepted_at': new['accepted_at'] = 'second'
    elif change == 'stable_field': new['target']['identity']['bytes'] += 1
    else: new['target']['identity']['inode'] += 1
    with pytest.raises(RuntimeError): a.audit_chain(old, new)
