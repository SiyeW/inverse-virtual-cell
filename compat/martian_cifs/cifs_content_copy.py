"""Explicit content-only copy for a report attachment on the approved output FS.

The destination inherits its filesystem permissions; source permissions are not
promised or simulated. This helper is opt-in, bounded and content-verified.
"""
import hashlib,os,shutil,stat
from pathlib import Path

def checked_copy(source,destination):
    if os.environ.get('MRO_CIFS_COMPAT')!='1':return shutil.copy(source,destination)
    root=Path(os.environ['MRO_CIFS_SAME_MODE_ROOT']).resolve(strict=True)
    src=Path(source);dst=Path(destination)
    if not dst.resolve().is_relative_to(root) or dst.is_symlink():raise ValueError('Report destination outside the owned output root')
    # Open first: CIFS may refresh a cached path timestamp on open. Compare the
    # same open file descriptor before/after reading, not a stale pre-open stat.
    with src.open('rb') as stream:
        before=os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size>20_000_000:raise ValueError('Not a bounded report attachment')
        data=stream.read(20_000_001);after=os.fstat(stream.fileno())
    h=hashlib.sha256(data).digest()
    if len(data)!=before.st_size or len(data)>20_000_000:raise ValueError('Report source size changed')
    if (before.st_size,before.st_mtime_ns,before.st_dev,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_dev,after.st_ino):
        raise ValueError('Report source identity changed')
    if dst.exists():
        if dst.stat().st_size!=len(data) or hashlib.sha256(dst.read_bytes()).digest()!=h:raise ValueError('Existing report destination differs')
    else:
        with dst.open('xb') as out:out.write(data);out.flush();os.fsync(out.fileno())
        if hashlib.sha256(dst.read_bytes()).digest()!=h:raise ValueError('Report attachment hash mismatch')
    return str(dst)
