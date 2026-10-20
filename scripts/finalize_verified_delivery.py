"""Finalize complete content receipts on CIFS without restarting a live copy.

This mode never copies or hashes a data file. Only a destination inode-only
remap can reuse an existing full SHA256/MD5 receipt; every other change stops.
Private manifests, handoff records and receipts must remain outside Git.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time

import copy_verified_delivery as delivery

POLICY = 'cifs-inode-audit-v1'
IDENTITY_FIELDS = {'bytes', 'mtime_ns', 'inode', 'device'}
STABLE_FIELDS = IDENTITY_FIELDS - {'inode'}


def read_manifest(path, digest, control):
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError('Frozen manifest SHA256 mismatch')
    manifest = json.loads(raw)
    delivery.validate_manifest(manifest)
    if Path(control).resolve() != Path(manifest['control']).resolve():
        raise ValueError('Control directory mismatch')
    if not set(manifest['vendor_md5']).issubset(manifest['files']):
        raise ValueError('Vendor hash outside frozen scope')
    return manifest


def require_cifs(target):
    result = subprocess.run(
        ['findmnt', '-n', '-T', str(target), '-o', 'FSTYPE'],
        capture_output=True, text=True, timeout=10, check=True)
    if result.stdout.strip() != 'cifs':
        raise RuntimeError('Inode-audit policy is restricted to a confirmed CIFS target')


def validate_identity(value):
    if set(value) != IDENTITY_FIELDS or any(type(v) is not int for v in value.values()):
        raise ValueError('Invalid receipt identity')


def validate_receipts(manifest, ledger, digest):
    """Validate proof values, not fresh file content, before any handoff signal."""
    if ledger.get('manifest_sha256') != digest:
        raise ValueError('Ledger belongs to another manifest')
    if set(ledger.get('verified', {})) != set(manifest['files']):
        raise RuntimeError('Finalize-only requires every frozen content receipt')
    for rel, frozen in manifest['files'].items():
        proof = ledger['verified'][rel]
        for side in ('source', 'target'):
            value = proof[side]
            validate_identity(value['identity'])
            if not re.fullmatch('[0-9a-f]{64}', value.get('sha256', '')):
                raise ValueError(f'Invalid SHA256 receipt: {rel}')
            if not re.fullmatch('[0-9a-f]{32}', value.get('md5', '')):
                raise ValueError(f'Invalid MD5 receipt: {rel}')
            if value['identity']['bytes'] != frozen['bytes']:
                raise RuntimeError(f'Receipt byte count changed: {rel}')
        if proof['source']['identity'] != frozen:
            raise RuntimeError(f'Frozen source receipt changed: {rel}')
        if any(proof['source'][k] != proof['target'][k] for k in ('sha256', 'md5')):
            raise RuntimeError(f'Content receipts disagree: {rel}')
        vendor = manifest['vendor_md5'].get(rel)
        if proof.get('vendor_md5') != vendor or (vendor and proof['source']['md5'] != vendor):
            raise RuntimeError(f'Vendor MD5 receipt changed: {rel}')
        if not proof.get('accepted_at'):
            raise ValueError(f'Missing completed publication receipt: {rel}')
        if 'target_content_identity' in proof:
            content_identity = proof['target_content_identity']
            validate_identity(content_identity)
            if any(content_identity[k] != proof['target']['identity'][k] for k in STABLE_FIELDS):
                raise RuntimeError(f'Audit changed content identity fields: {rel}')
        audits = proof.get('target_inode_audits', [])
        for audit in audits:
            validate_identity(audit['previous'])
            validate_identity(audit['current'])
            if (audit.get('policy') != POLICY or
                    any(audit['previous'][k] != audit['current'][k] for k in STABLE_FIELDS) or
                    any(audit.get(k) != proof['target'][k] for k in ('sha256', 'md5'))):
                raise RuntimeError(f'Invalid inode-only audit: {rel}')
        if audits and audits[-1]['current'] != proof['target']['identity']:
            raise RuntimeError(f'Inode audit chain does not match receipt: {rel}')


def reconcile_inode(proof, current):
    validate_identity(current)
    previous = proof['target']['identity']
    if any(current[k] != previous[k] for k in STABLE_FIELDS):
        raise RuntimeError('Target size/mtime/device changed; retained for explicit review')
    if current['inode'] == previous['inode']:
        return False
    # Preserve the identity tied to the last full content read. Updating the
    # observed inode below is explicitly metadata reconciliation, not rehashing.
    proof.setdefault('target_content_identity', copy.deepcopy(previous))
    proof.setdefault('target_inode_audits', []).append(dict(
        previous=copy.deepcopy(previous), current=copy.deepcopy(current),
        time=delivery.now(), policy=POLICY,
        method='inode-only metadata audit; existing full content hashes retained',
        sha256=proof['target']['sha256'], md5=proof['target']['md5']))
    proof['target']['identity'] = copy.deepcopy(current)
    return True


def directory_names(root):
    result = set()
    for base, dirs, _ in os.walk(root, followlinks=False):
        for name in dirs:
            path = Path(base)/name
            if path.is_symlink():
                raise ValueError('Directory symlink in delivery')
            result.add(path.relative_to(root).as_posix())
    return result


def preserve_control(control, folder, owner_raw=None):
    """Archive small, owned durable records after the previous writer exits."""
    folder.mkdir(exist_ok=False)
    records = {}
    for name in ('manifest.json', 'ledger.json', 'status.json', 'run.stdout.log',
                 'run.stderr.log', 'final_acceptance.json'):
        path = control/name
        if not path.exists():
            continue
        if path.stat().st_size > 4 * 1024**2:
            raise RuntimeError(f'Control archive exceeds bounded metadata budget: {name}')
        raw = path.read_bytes()
        with (folder/name).open('xb') as out:
            out.write(raw)
        records[name] = dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    if owner_raw is not None:
        with (folder/'owner.before.json').open('xb') as out:
            out.write(owner_raw)
        records['owner.before.json'] = dict(
            bytes=len(owner_raw), sha256=hashlib.sha256(owner_raw).hexdigest())
    delivery.atomic_json(folder/'snapshot.json', records)
    return records


def finalize(manifest_path, digest, control):
    """No fall-through to copy/resume/full-data hashes, even with missing files."""
    control = Path(control)
    manifest = read_manifest(manifest_path, digest, control)
    source, target = Path(manifest['source']), Path(manifest['target'])
    require_cifs(target)
    lock = control/'owner.lock'
    fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, 'w') as out:
        json.dump(dict(pid=os.getpid(), started_at=delivery.now(),
                       manifest_sha256=digest, finalize_only=True,
                       target_identity_policy=POLICY), out)

    def status(stage, **extra):
        delivery.atomic_json(control/'status.json', dict(
            stage=stage, time=delivery.now(), pid=os.getpid(), finalize_only=True,
            target_identity_policy=POLICY, total_files=manifest['total_files'],
            total_bytes=manifest['total_bytes'], **extra))

    try:
        status('final_acceptance')
        ledger_path = control/'ledger.json'
        original_raw = ledger_path.read_bytes()
        ledger = json.loads(original_raw)
        validate_receipts(manifest, ledger, digest)
        expected_dirs = set(manifest['directories'])
        if (delivery.tree_files(source) != manifest['files'] or
                directory_names(source) != expected_dirs):
            raise RuntimeError('Frozen source identity/scope changed')
        actual = delivery.tree_files(target)
        if set(actual) != set(manifest['files']):
            raise RuntimeError('Incomplete destination, extra files or active partials')
        if not directory_names(target).issubset(expected_dirs):
            raise RuntimeError('Unexpected destination directories')
        for rel in expected_dirs:
            delivery.safe_path(target, rel).mkdir(parents=True, exist_ok=True)
        if directory_names(target) != expected_dirs:
            raise RuntimeError('Destination directory coverage changed')
        for rel, current in actual.items():
            reconcile_inode(ledger['verified'][rel], current)
        # A second metadata view detects meaningful changes during the audit.
        # Inode remapping itself may repeat; it must not cause an unbounded loop.
        final_view = delivery.tree_files(target)
        if set(final_view) != set(manifest['files']):
            raise RuntimeError('Destination coverage changed during finalization')
        for rel, current in final_view.items():
            reconcile_inode(ledger['verified'][rel], current)
        if (delivery.tree_files(source) != manifest['files'] or
                directory_names(source) != expected_dirs):
            raise RuntimeError('Frozen source changed during finalization')
        validate_receipts(manifest, ledger, digest)
        if ledger_path.read_bytes() != original_raw:
            raise RuntimeError('Another writer changed ledger during finalization')
        history = control/'finalize_history'
        history.mkdir(exist_ok=True)
        snapshot = history/(hashlib.sha256(original_raw).hexdigest()+'_'+str(time.time_ns()))
        preserve_control(control, snapshot)
        ledger['target_identity_policy'] = POLICY
        delivery.atomic_json(ledger_path, ledger)
        audited = [rel for rel, proof in ledger['verified'].items()
                   if proof.get('target_inode_audits')]
        acceptance = dict(
            complete=True, completed_at=delivery.now(), source_retained=True,
            manifest_sha256=digest, total_files=manifest['total_files'],
            total_bytes=manifest['total_bytes'], vendor_md5_passed=len(manifest['vendor_md5']),
            source=str(source), target=str(target), finalize_only=True,
            target_identity_policy=POLICY, data_files_rehashed=0,
            inode_audited_files=len(audited),
            inode_audited_bytes=sum(manifest['files'][r]['bytes'] for r in audited),
            content_proof='existing full source/target SHA256 and MD5 receipts',
            prior_ledger_sha256=hashlib.sha256(original_raw).hexdigest(),
            ledger_sha256=hashlib.sha256(ledger_path.read_bytes()).hexdigest())
        delivery.atomic_json(control/'final_acceptance.json', acceptance)
        status('complete', verified_files=manifest['total_files'],
               verified_bytes=manifest['total_bytes'], source_retained=True)
        return acceptance
    except Exception as exc:
        status('failed', error=repr(exc))
        raise
    finally:
        if json.loads(lock.read_text())['pid'] != os.getpid():
            raise RuntimeError('Owner lock changed; not removing another owner')
        lock.unlink()


def process_view(pid):
    root = Path('/proc')/str(pid)
    try:
        before = (root/'stat').read_text().rsplit(')', 1)[1].split()
        command = (root/'cmdline').read_bytes().rstrip(b'\0').split(b'\0')
        uid = root.stat().st_uid
        after = (root/'stat').read_text().rsplit(')', 1)[1].split()
    except FileNotFoundError:
        return None
    if before[19] != after[19]:
        raise RuntimeError('Process identity changed during inspection')
    if after[0] == 'Z':
        return None
    return dict(startticks=int(after[19]), uid=uid,
                command=[part.decode() for part in command])


def require_old_owner(control, manifest_path, digest, pid, startticks, script, started_at):
    raw = (control/'owner.lock').read_bytes()
    owner = json.loads(raw)
    view = process_view(pid)
    expected_uid = os.getuid()
    if (owner.get('pid') != pid or owner.get('started_at') != started_at or
            owner.get('manifest_sha256') != digest or not view or
            view['startticks'] != startticks or view['uid'] != expected_uid or
            str(script) not in view['command']):
        raise RuntimeError('Previous owner PID/UID/start/command/lock changed')
    for flag, expected in (('--manifest', str(manifest_path)),
                           ('--manifest-sha256', digest), ('--control', str(control))):
        args = view['command']
        if args.count(flag) != 1 or args.index(flag)+1 >= len(args) or args[args.index(flag)+1] != expected:
            raise RuntimeError('Previous owner command scope changed')
    return raw


def handoff(manifest_path, digest, control, pid, startticks, script, started_at,
            operation, timeout=14400, stop_timeout=30, poll=1):
    """One-shot observer; signal only after every full content receipt exists."""
    control, operation = Path(control), Path(operation)
    manifest = read_manifest(manifest_path, digest, control)
    require_cifs(manifest['target'])
    if operation.parent.resolve() != control.resolve():
        raise ValueError('Handoff operation must be a new child of the control directory')
    operation.mkdir(exist_ok=False)
    record = dict(stage='waiting_for_all_content_receipts', started_at=delivery.now(),
                  observer_pid=os.getpid(), observer_uid=os.getuid(),
                  previous_pid=pid, previous_startticks=startticks,
                  previous_script=str(script), previous_started_at=started_at,
                  wait_timeout_seconds=timeout,
                  manifest_sha256=digest, signalled=False)
    delivery.atomic_json(operation/'status.json', record)
    deadline = time.monotonic()+timeout
    owner_raw = None
    try:
        while True:
            state = json.loads((control/'status.json').read_text())
            if state.get('stage') == 'failed':
                raise RuntimeError(f'Previous worker failed; no automatic recovery: {state.get("error")}')
            exists = (control/'owner.lock').exists()
            if exists:
                owner_raw = require_old_owner(control, manifest_path, digest, pid,
                                              startticks, script, started_at)
            elif (process_view(pid) is not None and
                  state.get('verified_files') != manifest['total_files']):
                raise RuntimeError('Previous owner lost its lock while still alive')
            if state.get('verified_files') == manifest['total_files']:
                ledger_raw = (control/'ledger.json').read_bytes()
                ledger = json.loads(ledger_raw)
                validate_receipts(manifest, ledger, digest)
                record.update(content_receipt_files=len(ledger['verified']),
                              ledger_before_signal_sha256=hashlib.sha256(ledger_raw).hexdigest())
                if exists and state.get('stage') != 'complete':
                    require_old_owner(control, manifest_path, digest, pid,
                                      startticks, script, started_at)
                    os.kill(pid, signal.SIGINT)
                    record['signalled'] = True
                    record['signal_reason'] = 'all content receipts durable; prevent redundant final rereads'
                break
            if not exists:
                raise RuntimeError('Previous writer exited before all content receipts; no restart')
            if time.monotonic() >= deadline:
                raise TimeoutError('Content gate deadline reached; previous writer left untouched')
            time.sleep(poll)
        stop_deadline = time.monotonic()+stop_timeout
        while process_view(pid) is not None or (control/'owner.lock').exists():
            if time.monotonic() >= stop_deadline:
                raise TimeoutError('Previous writer has not exited/released lock; not taking ownership')
            time.sleep(min(poll, .2))
        state = json.loads((control/'status.json').read_text())
        if state.get('stage') == 'failed':
            raise RuntimeError('Previous worker reported a real failure; inspect before finalizing')
        stderr = control/'run.stderr.log'
        if stderr.exists():
            if stderr.stat().st_size > 4 * 1024**2:
                raise RuntimeError('Previous stderr exceeds bounded review budget')
            error_text = stderr.read_text(errors='replace').strip()
            intentional = (record['signalled'] and error_text.endswith('KeyboardInterrupt') and
                           not re.search(r'(?m)^\w*(?:Error|Exception):', error_text))
            if error_text and not intentional:
                raise RuntimeError('Previous stderr is not solely an intentional interrupt; inspect first')
        record['snapshots'] = preserve_control(control, operation/'before', owner_raw)
        record.update(stage='old_writer_exited', completed_at=delivery.now(),
                      stderr_note='KeyboardInterrupt, when signalled here, is intentional handoff')
        delivery.atomic_json(operation/'status.json', record)
        return record
    except Exception as exc:
        record.update(stage='failed', error=repr(exc), time=delivery.now())
        delivery.atomic_json(operation/'status.json', record)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--handoff-pid', type=int)
    parser.add_argument('--handoff-startticks', type=int)
    parser.add_argument('--handoff-script', type=Path)
    parser.add_argument('--handoff-started-at')
    parser.add_argument('--handoff-operation', type=Path)
    parser.add_argument('--wait-timeout', type=float, default=14400)
    args = parser.parse_args()
    handoff_args = (args.handoff_pid, args.handoff_startticks, args.handoff_script,
                    args.handoff_started_at, args.handoff_operation)
    if any(value is not None for value in handoff_args):
        if not all(value is not None for value in handoff_args):
            parser.error('All handoff identity/operation arguments are required together')
        handoff(args.manifest, args.manifest_sha256, args.control, *handoff_args,
                timeout=args.wait_timeout)
    result = finalize(args.manifest, args.manifest_sha256, args.control)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
