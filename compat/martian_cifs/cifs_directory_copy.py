"""Content-verified directory publication in the explicitly owned CIFS tree.

Directory modes/timestamps are inherited from the destination. Native and
non-opted-in calls keep the standard library behavior. File publication still
uses Cell Ranger's original hardlink/copy callback.
"""
import hashlib
import os
from pathlib import Path
import shutil
import stat
import subprocess


def _cifs(root):
    if os.name != 'posix':
        return False
    result = subprocess.run(['/usr/bin/stat', '-f', '-c', '%t', str(root)],
                            capture_output=True, text=True, check=True, timeout=5)
    return result.stdout.strip().lower() in {'ff534d42', 'fe534d42'}


def _digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        before = os.fstat(stream.fileno())
        for data in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(data)
        after = os.fstat(stream.fileno())
    if (before.st_size, before.st_mtime_ns, before.st_dev) != (
            after.st_size, after.st_mtime_ns, after.st_dev):
        raise ValueError('Publication source changed while hashing')
    return h.digest()


def copytree(source, destination, *, copy_function, dirs_exist_ok):
    if os.environ.get('MRO_CIFS_COMPAT') != '1':
        return shutil.copytree(source, destination, copy_function=copy_function,
                               dirs_exist_ok=dirs_exist_ok)
    root = Path(os.environ['MRO_CIFS_SAME_MODE_ROOT']).resolve(strict=True)
    src, dst = Path(os.fsdecode(source)), Path(os.fsdecode(destination))
    if not _cifs(root):
        return shutil.copytree(source, destination, copy_function=copy_function,
                               dirs_exist_ok=dirs_exist_ok)

    def check(path):
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError('Directory publication escapes the owned output root')

    check(src)
    check(dst)
    if not src.is_dir() or dst.resolve() == src.resolve() or dst.resolve().is_relative_to(src.resolve()):
        raise ValueError('An independent publication directory is required')
    # Inventory before writing; this policy is limited to small result bundles.
    entries = []
    total = 0
    def inventory(parent):
        nonlocal total
        check(parent)
        for p in sorted(parent.iterdir()):
            check(p)
            mode = p.lstat().st_mode
            if not (stat.S_ISDIR(mode) or stat.S_ISREG(mode)):
                raise ValueError('Unsupported publication member')
            entries.append((p.relative_to(src), stat.S_ISDIR(mode)))
            if len(entries) > 10000:
                raise ValueError('Publication entry limit exceeded')
            if stat.S_ISDIR(mode):
                inventory(p)
            else:
                total += p.stat().st_size
                if total > 100_000_000:
                    raise ValueError('Publication byte limit exceeded')
    inventory(src)
    dst.mkdir(exist_ok=dirs_exist_ok)
    for relative, directory in entries:
        p, q = src / relative, dst / relative
        check(p)
        check(q)
        if directory:
            q.mkdir(exist_ok=dirs_exist_ok)
            continue
        expected = _digest(p)
        if q.exists():
            if not q.is_file() or _digest(q) != expected:
                raise ValueError('Existing publication file differs')
        else:
            copy_function(os.fsencode(p) if isinstance(source, bytes) else str(p),
                          os.fsencode(q) if isinstance(destination, bytes) else str(q))
        if _digest(q) != expected or _digest(p) != expected:
            raise ValueError('Published content differs or source changed')
    return destination
