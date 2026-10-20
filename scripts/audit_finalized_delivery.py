"""Read-only, independent audit of receipt-only CIFS delivery acceptance.

Only small control records are read. Data files are opened for refreshed fstat,
never read or hashed. Inode observations after acceptance are recorded, not
silently substituted into the accepted ledger or called fresh content hashes.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import time

POLICY = 'cifs-inode-audit-v1'
FIELDS = {'bytes', 'mtime_ns', 'inode', 'device'}
STABLE = FIELDS - {'inode'}
RECORD_LIMIT = 4 * 1024**2


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def record(path):
    path = Path(path)
    info = path.lstat()
    require(stat.S_ISREG(info.st_mode) and info.st_size <= RECORD_LIMIT,
            f'Not a bounded regular control record: {path}')
    raw = path.read_bytes()
    require(len(raw) <= RECORD_LIMIT, 'Control record grew past budget')
    return raw


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def regular_identity(path):
    info = Path(path).lstat()
    require(stat.S_ISREG(info.st_mode), f'Not a regular data file: {path}')
    flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0)
    fd = os.open(path, flags)
    try:
        info = os.fstat(fd)  # Refresh CIFS metadata; no data read.
        require(stat.S_ISREG(info.st_mode), 'Data file type changed')
        return dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns,
                    inode=info.st_ino, device=info.st_dev)
    finally:
        os.close(fd)


def inventory(root):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), 'Invalid inventory root')
    files, directories = {}, set()
    pending = [root]
    while pending:
        base = pending.pop()
        with os.scandir(base) as entries:
            for entry in entries:
                path = Path(entry.path)
                rel = path.relative_to(root).as_posix()
                require(not entry.is_symlink(), f'Symlink in delivery: {rel}')
                if entry.is_dir(follow_symlinks=False):
                    directories.add(rel)
                    pending.append(path)
                else:
                    files[rel] = regular_identity(path)
    return files, directories


def valid_identity(value):
    require(set(value) == FIELDS and all(type(v) is int for v in value.values()),
            'Invalid exact-integer file identity')


def process_identity(pid):
    root = Path('/proc') / str(pid)
    try:
        before = (root/'stat').read_text().rsplit(')', 1)[1].split()
        command = (root/'cmdline').read_bytes().replace(bytes([0]), b' ').decode()
        after = (root/'stat').read_text().rsplit(')', 1)[1].split()
    except FileNotFoundError:
        return None
    require(before[19] == after[19], 'PID changed during read-only inspection')
    return dict(startticks=int(after[19]), state=after[0], command=command,
                uid=root.stat().st_uid)


def mount_type(path):
    return subprocess.run(['findmnt', '-n', '-T', str(path), '-o', 'FSTYPE'],
                          check=True, capture_output=True, text=True,
                          timeout=10).stdout.strip()


def audit_chain(old, new):
    """Prove that only explicitly audited inode fields changed in each proof."""
    old_events = old.get('target_inode_audits', [])
    events = new.get('target_inode_audits', [])
    require(events[:len(old_events)] == old_events, 'Previous inode audit was changed')
    cursor = old['target']['identity']
    for event in events[len(old_events):]:
        valid_identity(event['previous'])
        valid_identity(event['current'])
        require(event['previous'] == cursor, 'Disconnected inode audit chain')
        require(event.get('policy') == POLICY and event.get('time') and
                event.get('method') == 'inode-only metadata audit; existing full content hashes retained',
                'Unlabelled metadata-only reconciliation')
        require(all(event['previous'][k] == event['current'][k] for k in STABLE),
                'Non-inode field changed inside audit')
        require(all(event.get(k) == old['target'][k] for k in ('sha256', 'md5')),
                'Audit changed original content hashes')
        cursor = event['current']
    require(cursor == new['target']['identity'], 'Audit does not reach accepted identity')
    normalized = copy.deepcopy(new)
    normalized['target']['identity'] = copy.deepcopy(old['target']['identity'])
    for key in ('target_inode_audits', 'target_content_identity'):
        if key in old:
            normalized[key] = copy.deepcopy(old[key])
        else:
            normalized.pop(key, None)
    require(normalized == old, 'Original content proof/publication fields were changed')
    if events:
        require(new.get('target_content_identity') == old.get(
            'target_content_identity', old['target']['identity']),
            'Original full-content identity was not preserved')


def audit(control, manifest_sha256, handoff_name, expected_exited=()):
    started = time.monotonic()
    control = Path(control)
    require(not (control/'owner.lock').exists(), 'Production owner still active')
    raw_manifest = record(control/'manifest.json')
    require(digest(raw_manifest) == manifest_sha256, 'Frozen manifest SHA256 changed')
    m = json.loads(raw_manifest)
    require(m.get('cleanup_sources') is False and Path(m['control']).resolve() == control.resolve(),
            'Wrong control scope/source-retention contract')
    require(len(m['files']) == m['total_files'] and
            sum(v['bytes'] for v in m['files'].values()) == m['total_bytes'],
            'Manifest totals disagree')
    for rel in list(m['files']) + list(m['directories']):
        p = PurePosixPath(rel)
        require(rel and not p.is_absolute() and '..' not in p.parts and '\\' not in rel,
                'Unsafe frozen relative path')
    require(set(m['vendor_md5']).issubset(m['files']), 'Vendor proof outside scope')
    source, target = Path(m['source']), Path(m['target'])
    require(mount_type(target) == 'cifs', 'Inode-only policy requires actual CIFS')
    originals = {name: record(control/name) for name in
                 ('ledger.json', 'final_acceptance.json', 'status.json')}
    ledger, acceptance, status_value = [json.loads(originals[n]) for n in
                                      ('ledger.json', 'final_acceptance.json', 'status.json')]
    require(status_value.get('stage') == 'complete' and status_value.get('finalize_only') is True,
            'Receipt-only finalization has not completed')
    require(acceptance.get('complete') is True and acceptance.get('finalize_only') is True and
            acceptance.get('source_retained') is True and acceptance.get('data_files_rehashed') == 0,
            'Missing receipt-only acceptance/source-retention proof')
    for k in ('total_files', 'total_bytes', 'source', 'target'):
        require(acceptance[k] == m[k], f'Acceptance scope differs: {k}')
    require(acceptance.get('vendor_md5_passed') == len(m['vendor_md5']) and
            acceptance.get('manifest_sha256') == manifest_sha256 and
            ledger.get('manifest_sha256') == manifest_sha256 and
            acceptance.get('target_identity_policy') == ledger.get('target_identity_policy') == POLICY and
            acceptance.get('ledger_sha256') == digest(originals['ledger.json']),
            'Acceptance/manifest/ledger binding failed')
    require(set(ledger['verified']) == set(m['files']), 'Incomplete/excess receipt coverage')

    require(handoff_name and Path(handoff_name).name == handoff_name, 'Unsafe handoff child')
    operation = control/handoff_name
    gate_raw = record(operation/'status.json')
    gate = json.loads(gate_raw)
    require(gate.get('stage') == 'old_writer_exited' and
            gate.get('manifest_sha256') == manifest_sha256 and
            gate.get('content_receipt_files') == m['total_files'], 'Handoff did not prove full content gate')
    snapshot_raw = record(operation/'before'/'snapshot.json')
    snapshots = json.loads(snapshot_raw)
    require('ledger.json' in snapshots and 'owner.before.json' in snapshots,
            'Handoff archive missing ledger/owner')
    for name, proof in snapshots.items():
        require(Path(name).name == name, 'Unsafe archive record name')
        raw = record(operation/'before'/name)
        require(proof == dict(bytes=len(raw), sha256=digest(raw)), 'Handoff archive SHA mismatch')
    prior_raw = record(operation/'before'/'ledger.json')
    prior_sha = digest(prior_raw)
    prior = json.loads(prior_raw)
    require(prior_sha == acceptance['prior_ledger_sha256'] == gate.get('ledger_before_signal_sha256'),
            'Pre-signal/prior-ledger archive binding changed')
    old_owner = json.loads(record(operation/'before'/'owner.before.json'))
    require(old_owner.get('pid') == gate['previous_pid'] and
            old_owner.get('started_at') == gate['previous_started_at'] and
            old_owner.get('manifest_sha256') == manifest_sha256,
            'Handoff owner archive disagrees')
    require(set(prior['verified']) == set(m['files']), 'Prior content ledger was incomplete')
    with os.scandir(control/'finalize_history') as entries:
        candidates = [Path(e.path) for e in entries if e.name.startswith(prior_sha+'_')]
    require(len(candidates) == 1 and record(candidates[0]/'ledger.json') == prior_raw,
            'Preserved finalizer history does not match prior ledger')
    previous_error = record(operation/'before'/'run.stderr.log').decode(errors='replace').strip() \
        if 'run.stderr.log' in snapshots else ''
    if previous_error:
        require(gate.get('signalled') is True and previous_error.endswith('KeyboardInterrupt') and
                not re.search(r'(?m)^\w*(?:Error|Exception):', previous_error),
                'Previous stderr is not solely the recorded intentional handoff')
    finalizer_error = record(control/'finalizer.stderr.log').decode(errors='replace').strip() \
        if (control/'finalizer.stderr.log').exists() else ''
    require(not finalizer_error, 'Finalizer stderr requires review')

    groups, inode_observations = {}, []
    audited_files = audited_bytes = 0
    for rel, frozen in m['files'].items():
        valid_identity(frozen)
        proof = ledger['verified'][rel]
        for side in ('source', 'target'):
            valid_identity(proof[side]['identity'])
            require(bool(re.fullmatch('[0-9a-f]{64}', proof[side].get('sha256', ''))) and
                    bool(re.fullmatch('[0-9a-f]{32}', proof[side].get('md5', ''))),
                    'Malformed full-content hash receipt')
        require(proof['source']['identity'] == frozen and proof['accepted_at'] and
                proof['target']['identity']['bytes'] == frozen['bytes'] and
                all(proof['source'][k] == proof['target'][k] for k in ('sha256', 'md5')),
                'Source/target content receipt mismatch')
        vendor = m['vendor_md5'].get(rel)
        require(proof.get('vendor_md5') == vendor and
                (vendor is None or proof['source']['md5'] == vendor), 'Vendor MD5 mismatch')
        audit_chain(prior['verified'][rel], proof)
        if proof.get('target_inode_audits'):
            audited_files += 1
            audited_bytes += frozen['bytes']
        group = groups.setdefault(rel.split('/')[0], dict(files=0, bytes=0, vendor_md5=0))
        group['files'] += 1
        group['bytes'] += frozen['bytes']
        group['vendor_md5'] += vendor is not None
    require((audited_files, audited_bytes) ==
            (acceptance['inode_audited_files'], acceptance['inode_audited_bytes']),
            'Inode-audit totals disagree')
    current_source, source_dirs = inventory(source)
    require(current_source == m['files'] and source_dirs == set(m['directories']),
            'Current source identity/file/directory coverage changed')
    current_target, target_dirs = inventory(target)
    require(set(current_target) == set(m['files']) and target_dirs == set(m['directories']),
            'Target exact coverage failed; missing/extra/partial entry')
    for rel, current in current_target.items():
        accepted = ledger['verified'][rel]['target']['identity']
        require(all(current[k] == accepted[k] for k in STABLE),
                f'Target stable identity changed after acceptance: {rel}')
        if current['inode'] != accepted['inode']:
            inode_observations.append(dict(file=rel, accepted=accepted, observed=current))
    actors = []
    for pid, ticks in expected_exited:
        view = process_identity(pid)
        require(not view or view['state'] == 'Z' or view['startticks'] != ticks,
                f'Owned writer/observer still alive: {pid}')
        actors.append(dict(pid=pid, startticks=ticks, owned_process_exited=True,
                           reused_pid_observation=view))
    require(not (control/'owner.lock').exists(), 'Writer acquired lock during audit')
    for name, raw in originals.items():
        require(record(control/name) == raw, 'Acceptance records changed during independent audit')
    require(record(operation/'status.json') == gate_raw, 'Handoff record changed during audit')
    return dict(complete=True, independently_audited_at=datetime.now().astimezone().isoformat(),
                elapsed_seconds=time.monotonic()-started, total_files=m['total_files'],
                total_bytes=m['total_bytes'], vendor_md5_passed=len(m['vendor_md5']),
                manifest_sha256=manifest_sha256, final_acceptance_sha256=digest(originals['final_acceptance.json']),
                ledger_sha256=digest(originals['ledger.json']), prior_ledger_sha256=prior_sha,
                handoff_status_sha256=digest(gate_raw), handoff_archive_snapshot_sha256=digest(snapshot_raw),
                policy=POLICY, content_proof='retained complete source/target SHA256 and MD5 receipts',
                data_files_rehashed=0, data_body_bytes_read=0, source_identity_exact=True,
                exact_file_and_directory_coverage=True, no_active_partials=True,
                no_owner_lock=True, actors=actors, groups=groups,
                accepted_inode_audited_files=audited_files, accepted_inode_audited_bytes=audited_bytes,
                post_acceptance_inode_only_observations=inode_observations,
                method='Independent receipt, archive-chain and refreshed-metadata audit; not a fresh data hash')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--handoff', required=True)
    parser.add_argument('--expected-exited', action='append', default=[], metavar='PID:STARTTICKS')
    args = parser.parse_args()
    exited = [tuple(map(int, value.split(':'))) for value in args.expected_exited]
    print(json.dumps(audit(args.control, args.manifest_sha256, args.handoff, exited), indent=2))


if __name__ == '__main__':
    main()
