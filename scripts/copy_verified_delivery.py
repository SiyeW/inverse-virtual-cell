"""Copy a frozen delivery without source deletion; publish only hash-verified files.

Private manifests contain paths and vendor hashes and must stay outside Git.
Writes are sequential; prefix reads may overlap with bounded buffers.
A detached host owns the CLI process.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import time

CHUNK = 4 * 1024**2


def now():
    return datetime.now().astimezone().isoformat()


def identity(path):
    s = Path(path).lstat()
    if not stat.S_ISREG(s.st_mode):
        raise ValueError(f'Not a regular file: {path}')
    # CIFS may return cached client mtime after close/fsync. Opening for read
    # revalidates metadata against the server; keep strict identity checks on
    # that refreshed view rather than ignoring mtime or trusting file size.
    with Path(path).open('rb') as f:
        s = os.fstat(f.fileno())
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino, device=s.st_dev)


def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(value, f, ensure_ascii=False, indent=2)
        f.flush()
        os.fsync(f.fileno())
    # An SMB reader/lease can briefly deny replacing even an owned JSON file.
    # Retry only this metadata operation, for a bounded 1.55 seconds total.
    # Persistent permission errors still stop with the data partial retained.
    for attempt in range(6):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            if attempt == 5:
                raise
            time.sleep(.05 * 2**attempt)


def safe_path(root, rel):
    p = PurePosixPath(rel)
    if p.is_absolute() or not p.parts or '..' in p.parts or '\\' in rel:
        raise ValueError('Unsafe relative path')
    result = Path(root).joinpath(*p.parts)
    if Path(root).resolve() not in result.resolve().parents:
        raise ValueError('Path escapes root')
    # Reject links even when they resolve within the permitted tree.
    for part in [result, *list(result.parents)[:len(p.parts)-1]]:
        if part.is_symlink():
            raise ValueError('Symlink in delivery path')
    return result


def tree_files(root):
    root = Path(root)
    found = {}
    for base, dirs, files in os.walk(root, followlinks=False):
        for name in dirs:
            if (Path(base)/name).is_symlink():
                raise ValueError('Directory symlink in tree')
        for name in files:
            p = Path(base)/name
            found[p.relative_to(root).as_posix()] = identity(p)
    return found


def validate_manifest(m):
    if m.get('cleanup_sources') is not False:
        raise ValueError('This runner never cleans sources')
    src, dst = Path(m['source']), Path(m['target'])
    if src.is_symlink() or dst.is_symlink():
        raise ValueError('Tree root cannot be a symlink')
    a, b = src.resolve(), dst.resolve()
    if a == b or a in b.parents or b in a.parents:
        raise ValueError('Overlapping source and target')
    for rel in m['files']:
        safe_path(src, rel)
        safe_path(dst, rel)
        if rel + '.receiving' in m['files']:
            raise ValueError('Partial filename collides with source')
    if sum(v['bytes'] for v in m['files'].values()) != m['total_bytes']:
        raise ValueError('Invalid byte total')
    if len(m['files']) != m['total_files']:
        raise ValueError('Invalid file total')


def hashes(path, callback=lambda *args: None):
    before = identity(path)
    sha, md5, done = hashlib.sha256(), hashlib.md5(), 0
    with Path(path).open('rb') as f:
        while chunk := f.read(CHUNK):
            sha.update(chunk)
            md5.update(chunk)
            done += len(chunk)
            callback(done, before['bytes'])
    if identity(path) != before:
        raise RuntimeError('File changed during hash')
    return dict(sha256=sha.hexdigest(), md5=md5.hexdigest(), identity=before)


def verify_prefix(reader, prior, length, sha, md5, callback, parallel=False):
    """Compare every prefix byte; optional bounded read-ahead, hashes stay ordered."""
    if not length:
        return

    def accept(block, previous, expected, done):
        if len(block) != expected or block != previous:
            raise RuntimeError('Partial prefix mismatch; preserved for inspection')
        sha.update(block)
        md5.update(block)
        callback('verify_resume_prefix', done, length)

    if not parallel:
        done = 0
        while done < length:
            expected = min(CHUNK, length-done)
            block, previous = reader.read(expected), prior.read(expected)
            done += expected
            accept(block, previous, expected, done)
        return

    # One pending read per endpoint. At most two blocks are prefetched while
    # the previous pair is checked/hashed; no seek, append or hash reordering.
    with ThreadPoolExecutor(max_workers=2) as pool:
        done, expected = 0, min(CHUNK, length)
        source_read = pool.submit(reader.read, expected)
        partial_read = pool.submit(prior.read, expected)
        while done < length:
            block, previous = source_read.result(), partial_read.result()
            next_done = done + expected
            next_expected = min(CHUNK, length-next_done)
            if next_expected:
                source_read = pool.submit(reader.read, next_expected)
                partial_read = pool.submit(prior.read, next_expected)
            accept(block, previous, expected, next_done)
            done, expected = next_done, next_expected


def copy_one(source, target, frozen, vendor_md5=None, callback=lambda *args: None,
             parallel_prefix=False):
    """Retain failed partials; never overwrite an existing published target."""
    source, target = Path(source), Path(target)
    if identity(source) != frozen:
        raise RuntimeError('Frozen source changed')
    if target.exists() or target.is_symlink():
        raise FileExistsError('Published destination already exists')
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.receiving')
    offset = identity(partial)['bytes'] if partial.exists() else 0
    if partial.is_symlink() or offset > frozen['bytes']:
        raise RuntimeError('Invalid partial')
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with source.open('rb') as reader:
        if offset:
            # Compare actual prefix bytes, not size/mtime, before resuming writes.
            with partial.open('rb') as prior:
                verify_prefix(reader, prior, offset, sha, md5, callback, parallel_prefix)
        with partial.open('ab' if partial.exists() else 'xb') as writer:
            done = offset
            while chunk := reader.read(CHUNK):
                writer.write(chunk)
                sha.update(chunk)
                md5.update(chunk)
                done += len(chunk)
                callback('copy', done, frozen['bytes'])
            writer.flush()
            os.fsync(writer.fileno())
    if identity(source) != frozen or done != frozen['bytes']:
        raise RuntimeError('Source changed or short read')
    original = dict(sha256=sha.hexdigest(), md5=md5.hexdigest(), identity=frozen)
    if vendor_md5 and original['md5'] != vendor_md5:
        raise RuntimeError('Vendor MD5 mismatch; partial retained, not published')
    verified = hashes(partial, lambda n, total: callback('hash_target', n, total))
    if verified['sha256'] != original['sha256']:
        raise RuntimeError('Destination SHA256 mismatch; partial retained')
    if identity(source) != frozen or target.exists() or target.is_symlink():
        raise RuntimeError('Source or destination changed before publish')
    os.rename(partial, target)
    verified['identity'] = identity(target)
    return dict(source=original, target=verified, vendor_md5=vendor_md5, accepted_at=now())


def revalidate_published(source, target, frozen, proof, callback=lambda *args: None):
    """A remapped destination inode requires fresh content proof, never a waiver."""
    if identity(source) != frozen:
        raise RuntimeError('Frozen source changed during receipt reconciliation')
    current = identity(target)
    if current == proof['target']['identity']:
        return False
    fresh = hashes(target, callback)
    if (fresh['sha256'] != proof['source']['sha256'] or
            fresh['md5'] != proof['source']['md5'] or
            (proof.get('vendor_md5') and fresh['md5'] != proof['vendor_md5'])):
        raise RuntimeError('Published destination content changed; retained without overwrite')
    if identity(source) != frozen:
        raise RuntimeError('Source changed during receipt reconciliation')
    proof.setdefault('identity_revalidations', []).append(dict(
        previous=proof['target']['identity'], current=fresh['identity'], time=now(),
        sha256=fresh['sha256'], method='fresh full destination hash, frozen source retained'))
    proof['target'] = fresh
    return True


def run(manifest_path, manifest_sha256, control, parallel_prefix=False):
    raw = Path(manifest_path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest_sha256:
        raise ValueError('Frozen manifest SHA256 mismatch')
    m = json.loads(raw)
    validate_manifest(m)
    control = Path(control)
    if control.resolve() != Path(m['control']).resolve():
        raise ValueError('Control directory mismatch')
    src, dst = Path(m['source']), Path(m['target'])
    lock = control/'owner.lock'
    fd = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, 'w') as f:
        json.dump(dict(pid=os.getpid(), started_at=now(), manifest_sha256=manifest_sha256,
                       parallel_prefix=parallel_prefix), f)
    ledger_path = control/'ledger.json'
    ledger = {'manifest_sha256':manifest_sha256, 'verified':{}, 'started_at':now()}
    started, last_status = time.monotonic(), 0.0

    def status(stage, force=False, **extra):
        nonlocal last_status
        t = time.monotonic()
        if not force and t-last_status < 3:
            return
        last_status = t
        records = ledger['verified']
        atomic_json(control/'status.json', dict(stage=stage, time=now(), pid=os.getpid(),
            parallel_prefix=parallel_prefix,
            total_files=m['total_files'], total_bytes=m['total_bytes'],
            verified_files=len(records), verified_bytes=sum(m['files'][r]['bytes'] for r in records),
            elapsed_seconds=t-started, **extra))

    def require_source():
        if not src.is_dir() or tree_files(src) != m['files']:
            raise RuntimeError('Source scope changed or disk unavailable')

    try:
        if ledger_path.exists():
            ledger = json.loads(ledger_path.read_text())
            if ledger['manifest_sha256'] != manifest_sha256:
                raise ValueError('Ledger belongs to another manifest')
        status('preflight', True)
        require_source()
        dst.mkdir(parents=True, exist_ok=True)
        allowed = set(m['files']) | {r+'.receiving' for r in m['files']}
        if not set(tree_files(dst)).issubset(allowed):
            raise RuntimeError('Unexpected destination files')
        # CIFS can remap inode numbers. Any changed receipt identity needs a
        # fresh full content hash; never silently ignore a metadata difference.
        for rel, proof in ledger['verified'].items():
            callback = lambda n,t: status('revalidate_target', file=rel, file_bytes_done=n, file_bytes_total=t)
            if revalidate_published(safe_path(src,rel), safe_path(dst,rel), m['files'][rel], proof, callback):
                atomic_json(ledger_path, ledger)
        atomic_json(ledger_path, ledger)
        # Cell Ranger small files first: make scientific inputs available early.
        order = sorted(m['files'], key=lambda r:(not r.startswith('cellrangerresult/'), m['files'][r]['bytes'], r))
        for rel in order:
            source, target = safe_path(src, rel), safe_path(dst, rel)
            if rel in ledger['verified']:
                continue
            frozen = m['files'][rel]
            vendor = m['vendor_md5'].get(rel)
            callback = lambda stage, n, total: status(stage, file=rel, file_bytes_done=n, file_bytes_total=total)
            status('start_file', True, file=rel, file_bytes_done=0, file_bytes_total=frozen['bytes'])
            if target.exists():
                a = hashes(source, lambda n,t: callback('hash_existing_source',n,t))
                b = hashes(target, lambda n,t: callback('hash_existing_target',n,t))
                if a['identity'] != frozen or a['sha256'] != b['sha256'] or (vendor and a['md5'] != vendor):
                    raise RuntimeError('Existing destination or source mismatch')
                proof = dict(source=a,target=b,vendor_md5=vendor,accepted_at=now())
            else:
                proof = copy_one(source,target,frozen,vendor,callback,parallel_prefix)
            ledger['verified'][rel] = proof
            atomic_json(ledger_path, ledger)
            status('file_verified', True, file=rel)
        for rel in m['directories']:
            safe_path(dst, rel).mkdir(parents=True, exist_ok=True)
        status('final_acceptance', True)
        require_source()
        actual = tree_files(dst)
        if set(actual) != set(m['files']) or set(ledger['verified']) != set(m['files']):
            raise RuntimeError('Incomplete destination or extra partials')
        for rel, proof in ledger['verified'].items():
            callback = lambda n,t: status('revalidate_target', file=rel, file_bytes_done=n, file_bytes_total=t)
            if revalidate_published(safe_path(src,rel), safe_path(dst,rel), m['files'][rel], proof, callback):
                atomic_json(ledger_path, ledger)
            if identity(safe_path(dst,rel)) != proof['target']['identity'] or proof['source']['sha256'] != proof['target']['sha256']:
                raise RuntimeError('Final destination proof changed')
        require_source()
        acceptance = dict(complete=True, completed_at=now(), source_retained=True,
            parallel_prefix=parallel_prefix,
            manifest_sha256=manifest_sha256,total_files=m['total_files'],total_bytes=m['total_bytes'],
            vendor_md5_passed=len(m['vendor_md5']),source=str(src),target=str(dst),
            ledger_sha256=hashlib.sha256(ledger_path.read_bytes()).hexdigest())
        atomic_json(control/'final_acceptance.json', acceptance)
        status('complete', True, source_retained=True)
    except Exception as exc:
        status('failed', True, error=repr(exc))
        raise
    finally:
        lock.unlink()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--manifest-sha256', required=True)
    p.add_argument('--control', type=Path, required=True)
    p.add_argument('--parallel-prefix', action='store_true',
                   help='Overlap bounded source/partial reads; still compare the entire prefix')
    a = p.parse_args()
    run(a.manifest, a.manifest_sha256, a.control, a.parallel_prefix)
