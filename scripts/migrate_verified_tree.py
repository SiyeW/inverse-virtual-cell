"""Resumable Linux tree migration with hashes, frozen scope and exact cleanup.

Configuration and journal remain private operational artifacts outside Git.
No source cleanup occurs until every tree has passed source/destination SHA256
and any supplied vendor MD5. Compatibility symlinks live on the source filesystem.
"""
from __future__ import annotations
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time


def atomic_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',encoding='utf8') as f:
        json.dump(value,f,ensure_ascii=False,indent=2)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp,path)


def identity(path):
    s=path.stat()
    return {'bytes':s.st_size,'mtime_ns':s.st_mtime_ns,'inode':s.st_ino}


def descendant(path, base):
    path,base=path.resolve(),base.resolve()
    if path==base or base not in path.parents:
        raise ValueError(f'Target outside permitted root: {path}')
    return path


def inventory(root):
    if root.is_symlink() or not root.is_dir():
        raise ValueError('Source must be an actual directory')
    result={}
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise ValueError(f'Symlink inside source tree: {path.name}')
        if path.is_file():
            result[str(path.relative_to(root))]=identity(path)
        elif not path.is_dir():
            raise ValueError('Special file in source tree')
    return result


def vendor_hashes(root):
    result={}
    for manifest in sorted(root.rglob('md5.txt')):
        for line in manifest.read_text(errors='strict').splitlines():
            match=re.fullmatch(r'([0-9a-fA-F]{32})\s+\*?(.+)',line.strip())
            if not match:
                if line.strip():
                    raise ValueError('Malformed vendor MD5 line')
                continue
            path=descendant(manifest.parent/match[2],root)
            rel=str(path.relative_to(root.resolve()))
            if rel in result and result[rel]!=match[1].lower():
                raise ValueError('Conflicting MD5 records')
            result[rel]=match[1].lower()
    return result


def hash_file(path, callback=None):
    before=identity(path)
    sha,md5=hashlib.sha256(),hashlib.md5()
    count,last=0,time.monotonic()
    with path.open('rb') as f:
        while chunk:=f.read(4*1024**2):
            sha.update(chunk)
            md5.update(chunk)
            count+=len(chunk)
            if callback and time.monotonic()-last>=3:
                callback(count,before['bytes'])
                last=time.monotonic()
    if identity(path)!=before:
        raise RuntimeError('File changed during hashing')
    if callback:
        callback(count,before['bytes'])
    return {'sha256':sha.hexdigest(),'md5':md5.hexdigest(),'identity':before}


def plan(config):
    source_base=Path(config['source_base']).resolve()
    target_base=Path(config['target_base']).resolve()
    if source_base==target_base or source_base in target_base.parents or target_base in source_base.parents:
        raise ValueError('Source and destination roots overlap')
    sources,targets=[],[]
    groups=[]
    for entry in config['groups']:
        source=descendant(Path(entry['source']),source_base)
        target=descendant(target_base/entry['relative_target'],target_base)
        for old in sources:
            if source==old or source in old.parents or old in source.parents:
                raise ValueError('Source groups overlap')
        for old in targets:
            if target==old or target in old.parents or old in target.parents:
                raise ValueError('Destination groups overlap')
        sources.append(source)
        targets.append(target)
        files=inventory(source)
        expected=vendor_hashes(source)
        expected.update({k:v.lower() for k,v in entry.get('expected_md5',{}).items()})
        if not set(expected).issubset(files):
            raise ValueError('Expected vendor file is missing')
        if target.exists() and any(target.iterdir()):
            raise FileExistsError('Initial destination must be empty')
        groups.append({'name':entry['name'],'source':str(source),'target':str(target),
                       'files':files,'expected_md5':expected,'verified':{},'source_hashes':{},
                       'cleanup_complete':False})
    return groups


def require_source_scope(group, root=None, allow_missing=False):
    root=Path(root or group['source'])
    actual=inventory(root)
    expected=group['files']
    if (not allow_missing and set(actual)!=set(expected)) or not set(actual).issubset(expected):
        raise RuntimeError('Source inventory changed')
    for rel,ident in actual.items():
        if ident!=expected[rel]:
            raise RuntimeError(f'Source identity changed: {rel}')


def require_verified_target(group):
    target=Path(group['target'])
    actual=inventory(target)
    if set(actual)!=set(group['files']) or set(group['verified'])!=set(actual):
        raise RuntimeError('Destination coverage incomplete or contains extra files')
    for rel,ident in actual.items():
        verified=group['verified'][rel]
        original=group['source_hashes'][rel]
        if ident!=verified['identity'] or verified['sha256']!=original['sha256']:
            raise RuntimeError(f'Destination verification identity changed: {rel}')


def cleanup_group(group, source_base, save):
    """Delete only frozen files after fresh scope and destination checks."""
    require_verified_target(group)
    original=Path(group['source'])
    quarantine=Path(group.get('quarantine',str(original)+'.verified_migration_20261015'))
    descendant(quarantine,source_base)
    if group['cleanup_complete']:
        if not original.is_symlink() or original.resolve()!=Path(group['target']).resolve():
            raise RuntimeError('Completed compatibility link changed')
        return
    if not group.get('cleanup_started'):
        require_source_scope(group)
        if quarantine.exists():
            raise FileExistsError('Unexpected quarantine path')
        group.update(cleanup_started=True,quarantine=str(quarantine))
        save()
    if not quarantine.exists():
        if original.is_dir() and not original.is_symlink():
            require_source_scope(group)
            original.rename(quarantine)
        elif not original.is_symlink():
            raise RuntimeError('Missing source/quarantine during cleanup')
    if not original.exists() and not original.is_symlink():
        original.symlink_to(Path(group['target']),target_is_directory=True)
    if not original.is_symlink() or original.resolve()!=Path(group['target']).resolve():
        raise RuntimeError('Unexpected source compatibility path')
    if quarantine.exists():
        require_source_scope(group,quarantine,allow_missing=True)
        for rel in group['files']:
            path=quarantine/rel
            if path.exists():
                if path.is_symlink() or identity(path)!=group['files'][rel]:
                    raise RuntimeError('Source changed immediately before unlink')
                path.unlink()
        for directory in sorted((p for p in quarantine.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
            directory.rmdir()
        quarantine.rmdir()
    group['cleanup_complete']=True
    save()


def run(config_path, control):
    if os.name!='posix':
        raise RuntimeError('Execution supported only on Linux')
    config_bytes=config_path.read_bytes()
    config=json.loads(config_bytes)
    control.mkdir(parents=True,exist_ok=True)
    lock=control/'owner.lock'
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    os.write(fd,str(os.getpid()).encode())
    os.close(fd)
    ledger_path,status_path=control/'ledger.json',control/'status.json'
    journal=None
    def save():
        atomic_json(ledger_path,journal)
    def status(stage,**detail):
        atomic_json(status_path,{'stage':stage,'pid':os.getpid(),
                    'time':datetime.now().astimezone().isoformat(),**detail})
    try:
        status('preflight')
        digest=hashlib.sha256(config_bytes).hexdigest()
        if ledger_path.exists():
            journal=json.loads(ledger_path.read_text())
            if journal['config_sha256']!=digest:
                raise ValueError('Configuration changed after migration started')
        else:
            groups=plan(config)
            total=sum(i['bytes'] for g in groups for i in g['files'].values())
            Path(config['target_base']).mkdir(parents=True,exist_ok=True)
            if shutil.disk_usage(config['target_base']).free<total+1024**3:
                raise OSError('Destination free space insufficient')
            journal={'config_sha256':digest,'groups':groups,'total_bytes':total,
                     'total_files':sum(len(g['files']) for g in groups),
                     'all_verified':False,'complete':False}
            save()
        if not journal['all_verified']:
            for group in journal['groups']:
                source,target=Path(group['source']),Path(group['target'])
                require_source_scope(group)
                target.mkdir(parents=True,exist_ok=True)
                for rel in group['files']:
                    cached=group['source_hashes'].get(rel)
                    if not cached or cached['identity']!=identity(source/rel):
                        callback=lambda done,total:status('hash_source',group=group['name'],file=rel,file_bytes_done=done,file_bytes_total=total)
                        cached=hash_file(source/rel,callback)
                        expected=group['expected_md5'].get(rel)
                        if expected and cached['md5']!=expected:
                            raise RuntimeError(f'Vendor MD5 mismatch: {rel}')
                        group['source_hashes'][rel]=cached
                        save()
                status('copy',group=group['name'])
                with (control/(group['name']+'.rsync.log')).open('ab') as log:
                    subprocess.run(['rsync','-rt','--partial','--append-verify','--info=progress2',
                                    str(source)+'/',str(target)+'/'],stdout=log,stderr=subprocess.STDOUT,check=True)
                for rel in group['files']:
                    cached=group['verified'].get(rel)
                    if not cached or cached['identity']!=identity(target/rel):
                        callback=lambda done,total:status('hash_target',group=group['name'],file=rel,file_bytes_done=done,file_bytes_total=total)
                        cached=hash_file(target/rel,callback)
                        if cached['sha256']!=group['source_hashes'][rel]['sha256']:
                            raise RuntimeError(f'Source/destination SHA256 mismatch: {rel}')
                        group['verified'][rel]=cached
                        save()
                require_verified_target(group)
                require_source_scope(group)
            journal['all_verified']=True
            save()
        # Freshly recheck all four scopes before ANY source file is removed.
        for group in journal['groups']:
            require_verified_target(group)
            if not group.get('cleanup_started'):
                require_source_scope(group)
        status('verified_all_sources_retained',files=journal['total_files'],bytes=journal['total_bytes'])
        if config.get('cleanup_sources'):
            for group in journal['groups']:
                status('cleanup_verified_source',group=group['name'])
                cleanup_group(group,Path(config['source_base']),save)
        journal['complete']=True
        journal['completed_at']=datetime.now().astimezone().isoformat()
        save()
        root=Path(config['target_base'])
        atomic_json(root/'provenance'/'migration_20261015.json',journal)
        lines=['# Perturb-seq data inventory','',f"Migration completed: {journal['completed_at']}",'',
               'All files matched source/destination SHA256 and configured vendor MD5.',
               'Local Windows datasets were not moved. Legacy data3 paths are compatibility symlinks.','']
        for group in journal['groups']:
            lines.extend([f"## {group['name']}",'',f"Path: `{group['target']}`",f"Files: {len(group['files'])}",
                          f"Bytes: {sum(x['bytes'] for x in group['files'].values())}",''])
        (root/'DATA_INVENTORY.md').write_text('\n'.join(lines),encoding='utf8')
        status('complete',files=journal['total_files'],bytes=journal['total_bytes'],cleanup=bool(config.get('cleanup_sources')),
               destination_free_bytes=shutil.disk_usage(root).free)
    except Exception as exc:
        status('failed',error=repr(exc))
        raise
    finally:
        lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True)
    a=p.parse_args()
    run(a.config,a.control)
