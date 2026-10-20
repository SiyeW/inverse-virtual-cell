import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

SCRIPTS = Path(__file__).parents[1]/'scripts'
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location('finalizer', SCRIPTS/'finalize_verified_delivery.py')
f = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f)
d = f.delivery


def completed_delivery(tmp_path, monkeypatch):
    source, target, control = [tmp_path/name for name in ('source', 'target', 'control')]
    source.mkdir(); control.mkdir()
    (source/'cellrangerresult').mkdir()
    (source/'cellrangerresult/a').write_bytes(b'abc'*100)
    vendor = hashlib.md5(b'abc'*100).hexdigest()
    manifest = dict(source=str(source), target=str(target), control=str(control),
                    cleanup_sources=False, files=d.tree_files(source),
                    directories=['cellrangerresult'], total_files=1, total_bytes=300,
                    vendor_md5={'cellrangerresult/a': vendor})
    path = control/'manifest.json'; path.write_text(json.dumps(manifest))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    d.run(path, digest, control)
    monkeypatch.setattr(f, 'require_cifs', lambda _: None)
    return manifest, path, digest, control


def load_ledger(control):
    return json.loads((control/'ledger.json').read_text())


def save_ledger(control, ledger):
    d.atomic_json(control/'ledger.json', ledger)


def forbid_data_work(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('Finalize-only must not copy, resume or hash any data file')
    monkeypatch.setattr(d, 'hashes', forbidden)
    monkeypatch.setattr(d, 'copy_one', forbidden)


def test_inode_audit_preserves_content_proof_without_any_data_rehash(tmp_path, monkeypatch):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    ledger = load_ledger(control)
    proof = ledger['verified']['cellrangerresult/a']
    proof['target']['identity']['inode'] += 10
    original = copy.deepcopy(proof)
    save_ledger(control, ledger)
    forbid_data_work(monkeypatch)
    result = f.finalize(path, digest, control)
    actual = load_ledger(control)['verified']['cellrangerresult/a']
    assert result['complete'] and result['data_files_rehashed'] == 0
    assert result['inode_audited_files'] == 1 and result['inode_audited_bytes'] == 300
    assert actual['source'] == original['source'] and actual['accepted_at'] == original['accepted_at']
    assert actual['target_content_identity'] == original['target']['identity']
    assert actual['target']['sha256'] == original['target']['sha256']
    assert len(actual['target_inode_audits']) == 1
    assert not (control/'owner.lock').exists()
    assert d.tree_files(m['source']) == m['files']
    assert len(list((control/'finalize_history').iterdir())) == 1
    # Repeating metadata reconciliation does not duplicate audit events or fail
    # because the preserved previous ledger has already been archived.
    f.finalize(path, digest, control)
    assert len(load_ledger(control)['verified']['cellrangerresult/a']['target_inode_audits']) == 1
    assert len(list((control/'finalize_history').iterdir())) == 2


@pytest.mark.parametrize('field', ['bytes', 'mtime_ns', 'device'])
def test_other_target_changes_stop_without_rehash(field):
    previous = dict(bytes=3, mtime_ns=42, inode=2, device=1)
    proof = dict(target=dict(identity=previous, sha256='a'*64, md5='b'*32))
    current = dict(previous); current[field] += 1; current['inode'] += 10
    with pytest.raises(RuntimeError, match='size/mtime/device'):
        f.reconcile_inode(proof, current)
    assert proof == dict(target=dict(identity=previous, sha256='a'*64, md5='b'*32))


def test_frozen_source_change_stops_and_leaves_ledger_untouched(tmp_path, monkeypatch):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    before = (control/'ledger.json').read_bytes()
    source = Path(m['source'])/'cellrangerresult/a'
    os.utime(source, ns=(source.stat().st_atime_ns, source.stat().st_mtime_ns+10**9))
    forbid_data_work(monkeypatch)
    with pytest.raises(RuntimeError, match='Frozen source'):
        f.finalize(path, digest, control)
    assert (control/'ledger.json').read_bytes() == before
    assert not (control/'owner.lock').exists()


@pytest.mark.parametrize('broken', ['sha256', 'md5', 'vendor', 'source_identity', 'accepted_at', 'missing'])
def test_bad_or_missing_receipts_never_fall_back_to_copy(tmp_path, monkeypatch, broken):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    ledger = load_ledger(control); proof = ledger['verified']['cellrangerresult/a']
    if broken in ('sha256', 'md5'):
        proof['target'][broken] = '0'*(64 if broken == 'sha256' else 32)
    elif broken == 'vendor': proof['vendor_md5'] = '0'*32
    elif broken == 'source_identity': proof['source']['identity']['inode'] += 1
    elif broken == 'accepted_at': proof.pop('accepted_at')
    else: ledger['verified'].clear()
    save_ledger(control, ledger); before = (control/'ledger.json').read_bytes()
    forbid_data_work(monkeypatch)
    with pytest.raises((RuntimeError, ValueError)):
        f.finalize(path, digest, control)
    assert (control/'ledger.json').read_bytes() == before
    assert (Path(m['target'])/'cellrangerresult/a').read_bytes() == b'abc'*100


@pytest.mark.parametrize('extra', ['partial', 'extra_dir', 'source_dir', 'missing_target'])
def test_exact_coverage_required(tmp_path, monkeypatch, extra):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    target = Path(m['target'])
    if extra == 'partial': (target/'extra.receiving').write_bytes(b'retain')
    elif extra == 'extra_dir': (target/'extra').mkdir()
    elif extra == 'source_dir': (Path(m['source'])/'extra').mkdir()
    else: (target/'cellrangerresult/a').unlink()
    forbid_data_work(monkeypatch)
    with pytest.raises(RuntimeError): f.finalize(path, digest, control)
    assert not (control/'owner.lock').exists()


def test_non_cifs_target_and_second_owner_rejected(tmp_path, monkeypatch):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    (control/'owner.lock').write_text('existing')
    with pytest.raises(FileExistsError): f.finalize(path, digest, control)
    assert (control/'owner.lock').read_text() == 'existing'
    monkeypatch.setattr(f.subprocess, 'run', lambda *a, **k: type('R', (), {'stdout': 'ext4\n'})())
    # Call the original function, not the test fixture's no-op mount checker.
    other_spec = importlib.util.spec_from_file_location('mountcheck', SCRIPTS/'finalize_verified_delivery.py')
    other = importlib.util.module_from_spec(other_spec); other_spec.loader.exec_module(other)
    with pytest.raises(RuntimeError, match='CIFS'): other.require_cifs(m['target'])


def test_second_metadata_inode_change_is_audited_without_a_retry_loop(tmp_path, monkeypatch):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    original = d.tree_files; calls = []
    def remapping(root):
        result = original(root)
        if Path(root) == Path(m['target']):
            calls.append(1)
            for identity in result.values(): identity['inode'] += len(calls)
        return result
    monkeypatch.setattr(d, 'tree_files', remapping)
    forbid_data_work(monkeypatch)
    result = f.finalize(path, digest, control)
    assert len(calls) == 2 and result['inode_audited_files'] == 1
    assert len(load_ledger(control)['verified']['cellrangerresult/a']['target_inode_audits']) == 2


def fake_old_owner(control, path, digest, monkeypatch):
    pid, ticks, started = 12345, 42, 'original-start'
    script = Path('/owned/copy_verified_delivery.py')
    (control/'owner.lock').write_text(json.dumps(dict(pid=pid, started_at=started, manifest_sha256=digest)))
    d.atomic_json(control/'status.json', dict(stage='final_acceptance', verified_files=1))
    view = dict(uid=2601, startticks=ticks, command=['python', str(script),
                '--manifest', str(path), '--manifest-sha256', digest, '--control', str(control)])
    monkeypatch.setattr(f.os, 'getuid', lambda: 2601, raising=False)
    state = {'alive': True, 'signals': []}
    monkeypatch.setattr(f, 'process_view', lambda _: view if state['alive'] else None)
    def interrupt(given_pid, sig):
        assert given_pid == pid and sig == signal.SIGINT
        assert len(load_ledger(control)['verified']) == 1
        state['signals'].append(sig); state['alive'] = False
        (control/'owner.lock').unlink()
        (control/'run.stderr.log').write_text('Traceback (most recent call last):\nKeyboardInterrupt\n')
    monkeypatch.setattr(f.os, 'kill', interrupt)
    return pid, ticks, script, started, view, state


def test_handoff_waits_for_durable_full_receipts_then_releases_writer(tmp_path, monkeypatch):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    pid, ticks, script, started, view, state = fake_old_owner(control, path, digest, monkeypatch)
    d.atomic_json(control/'status.json', dict(stage='hash_target', verified_files=0))
    def finish(_): d.atomic_json(control/'status.json', dict(stage='file_verified', verified_files=1))
    monkeypatch.setattr(f.time, 'sleep', finish)
    result = f.handoff(path, digest, control, pid, ticks, script, started, control/'handoff', timeout=3)
    assert state['signals'] == [signal.SIGINT] and result['stage'] == 'old_writer_exited'
    assert result['snapshots']['ledger.json']['sha256'] == hashlib.sha256((control/'ledger.json').read_bytes()).hexdigest()
    assert not (control/'owner.lock').exists()
    forbid_data_work(monkeypatch)
    assert f.finalize(path, digest, control)['complete']


@pytest.mark.parametrize('change', ['uid', 'startticks', 'script', 'scope', 'lock'])
def test_handoff_never_signals_wrong_or_reused_owner(tmp_path, monkeypatch, change):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    pid, ticks, script, started, view, state = fake_old_owner(control, path, digest, monkeypatch)
    if change in ('uid', 'startticks'): view[change] += 1
    elif change == 'script': view['command'][1] = '/unrelated/program.py'
    elif change == 'scope': view['command'][-1] = '/different/control'
    else: (control/'owner.lock').write_text(json.dumps(dict(pid=999, started_at=started, manifest_sha256=digest)))
    with pytest.raises(RuntimeError):
        f.handoff(path, digest, control, pid, ticks, script, started, control/'handoff')
    assert state['signals'] == [] and (control/'owner.lock').exists()


@pytest.mark.parametrize('failure', ['timeout', 'worker_failed', 'bad_receipt', 'unexplained_stderr'])
def test_handoff_failure_preserves_data_and_does_not_finalize(tmp_path, monkeypatch, failure):
    m, path, digest, control = completed_delivery(tmp_path, monkeypatch)
    pid, ticks, script, started, view, state = fake_old_owner(control, path, digest, monkeypatch)
    if failure == 'timeout': d.atomic_json(control/'status.json', dict(stage='copy', verified_files=0))
    elif failure == 'worker_failed': d.atomic_json(control/'status.json', dict(stage='failed', verified_files=1, error='content mismatch'))
    elif failure == 'bad_receipt':
        ledger = load_ledger(control); ledger['verified'].clear(); save_ledger(control, ledger)
    else:
        def interrupted_with_error(*args):
            state['signals'].append(signal.SIGINT); state['alive'] = False
            (control/'owner.lock').unlink()
            (control/'run.stderr.log').write_text('RuntimeError: actual content mismatch\nKeyboardInterrupt\n')
        monkeypatch.setattr(f.os, 'kill', interrupted_with_error)
    with pytest.raises((RuntimeError, TimeoutError)):
        f.handoff(path, digest, control, pid, ticks, script, started, control/'handoff', timeout=0)
    if failure != 'unexplained_stderr': assert state['signals'] == []
    assert (Path(m['target'])/'cellrangerresult/a').read_bytes() == b'abc'*100
    assert json.loads((control/'handoff/status.json').read_text())['stage'] == 'failed'


@pytest.mark.skipif(not os.environ.get('DELIVERY_CIFS_PROBE_ROOT'), reason='explicit real-CIFS probe root required')
def test_real_cifs_owned_process_handoff_and_receipt_only_finalization(tmp_path, monkeypatch):
    root = Path(os.environ['DELIVERY_CIFS_PROBE_ROOT'])
    root.mkdir(exist_ok=False)
    source, target, control = tmp_path/'source', root/'target', root/'control'
    source.mkdir(); control.mkdir(); (source/'cellrangerresult').mkdir()
    data = b'actual-cifs-content-proof'*1024
    (source/'cellrangerresult/a').write_bytes(data)
    manifest = dict(source=str(source), target=str(target), control=str(control),
                    cleanup_sources=False, files=d.tree_files(source), directories=['cellrangerresult'],
                    total_files=1, total_bytes=len(data), vendor_md5={'cellrangerresult/a': hashlib.md5(data).hexdigest()})
    path = control/'manifest.json'; path.write_text(json.dumps(manifest))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    worker = tmp_path/'owned_test_worker.py'
    worker.write_text('''import sys, json, os, time
from pathlib import Path
sys.path.insert(0, sys.argv[1])
import copy_verified_delivery as d
p=Path(sys.argv[sys.argv.index('--manifest')+1]); m=json.loads(p.read_text())
c=Path(m['control']); h=sys.argv[sys.argv.index('--manifest-sha256')+1]
(c/'owner.lock').write_text(json.dumps(dict(pid=os.getpid(),started_at='probe',manifest_sha256=h)))
try:
 proof=d.copy_one(Path(m['source'])/'cellrangerresult/a',Path(m['target'])/'cellrangerresult/a',m['files']['cellrangerresult/a'],m['vendor_md5']['cellrangerresult/a'])
 proof['target']['identity']['inode']+=1
 d.atomic_json(c/'ledger.json',dict(manifest_sha256=h,verified={'cellrangerresult/a':proof}))
 d.atomic_json(c/'status.json',dict(stage='final_acceptance',verified_files=1))
 while True: time.sleep(.1)
finally: (c/'owner.lock').unlink()
''')
    stderr = (control/'run.stderr.log').open('w')
    process = subprocess.Popen([sys.executable, str(worker), str(SCRIPTS), '--manifest', str(path),
                                '--manifest-sha256', digest, '--control', str(control)], stderr=stderr)
    try:
        deadline = time.monotonic()+60
        while not (control/'status.json').exists():
            if time.monotonic() >= deadline or process.poll() is not None:
                pytest.fail('Owned CIFS probe did not reach its content gate')
            time.sleep(.1)
        view = f.process_view(process.pid)
        result = f.handoff(path, digest, control, process.pid, view['startticks'], worker,
                           'probe', control/'handoff', timeout=10)
        assert result['signalled'] and process.wait(timeout=5) != 0
        forbid_data_work(monkeypatch)
        final = f.finalize(path, digest, control)
        assert final['complete'] and final['data_files_rehashed'] == 0 and final['inode_audited_files'] == 1
        assert (target/'cellrangerresult/a').read_bytes() == data
        assert (source/'cellrangerresult/a').read_bytes() == data
        assert not (control/'owner.lock').exists()
    finally:
        if process.poll() is None:
            process.terminate(); process.wait(timeout=5)
        stderr.close()
