"""Copy an explicit stopped-checkpoint inventory; preserve source and failed payloads.

No pipeline execution. A separate tested launch gate must admit the result.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
from clone_metadata_for_inspect import clone

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()

def identity(p):
    s=p.stat();return dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)

def atomic(p,obj):
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n')
    for i in range(6):
        try:tmp.replace(p);return
        except PermissionError:
            if i==5:raise
            time.sleep(.05*2**i)

def mappings(receipt):
    src=receipt['source'];dst=receipt['target']
    aliases=[(src+'/'+x['source_target'],dst+'/'+x['path']) for x in receipt['aliases_flattened']]
    return sorted(aliases,key=lambda x:len(x[0]),reverse=True)+[(src,dst)]

def mapped_path(p,maps):
    for old,new in maps:
        if p==old or p.startswith(old+'/'):return new+p[len(old):]
    raise ValueError('Path is outside the stopped source')

def rebase(data,maps):
    for old,new in maps:data=data.replace(old.encode(),new.encode())
    return data

def migrate(plan,control):
    source=Path(plan['source']).resolve(strict=True);target=Path(plan['target']).absolute()
    if str(source)!=plan['source'] or target.exists() or control.exists():raise ValueError('Fresh independent targets and canonical source required')
    if control.resolve().is_relative_to(source) or source.is_relative_to(control.resolve()):raise ValueError('Control overlaps source')
    files=plan['files']
    if len(files)!=plan['expected_files'] or sum(x['identity']['bytes'] for x in files)!=plan['expected_bytes']:raise ValueError('Inventory totals differ')
    if len({x['path'] for x in files})!=len(files):raise ValueError('Duplicate inventory path')
    for item in files:
        p=Path(item['path'])
        if not p.resolve(strict=True).is_relative_to(source) or p.is_symlink() or not p.is_file():raise ValueError('Non-owned regular payload')
        if identity(p)!=item['identity']:raise ValueError('Source identity differs: '+str(p))
    target.parent.mkdir(parents=True,exist_ok=True);control.mkdir(parents=True)
    rec=dict(stage='copying_metadata',pid=os.getpid(),uid=os.getuid() if hasattr(os,'getuid') else None,
             source=str(source),target=str(target),expected_files=len(files),expected_bytes=plan['expected_bytes'],
             accepted_files=[],copied_bytes=0,science_started=False)
    start=time.monotonic()
    def save():
        rec['elapsed_seconds']=time.monotonic()-start;atomic(control/'status.json',rec)
    with (control/'owner.lock').open('x') as f:json.dump(dict(pid=os.getpid()),f)
    atomic(control/'plan.json',plan);save()
    try:
        meta=clone(source,target);maps=mappings(meta)
        destinations=[mapped_path(x['path'],maps) for x in files]
        if len(set(destinations))!=len(files):raise ValueError('Collapsed destination conflict')
        for item,dest in zip(files,destinations):
            src=Path(item['path']);dst=Path(dest)
            if not dst.is_relative_to(target):raise ValueError('Escaping destination')
            dst.parent.mkdir(parents=True,exist_ok=True)
            part=dst.with_name(dst.name+'.checkpoint-copying')
            if dst.exists() or part.exists():raise ValueError('Payload destination occupied')
            rec.update(stage='copying_payload',current_file=str(src),file_bytes_done=0,file_bytes=item['identity']['bytes']);save()
            h=hashlib.sha256();last=time.monotonic()
            with src.open('rb') as f,part.open('xb') as g:
                if identity(src)!=item['identity']:raise ValueError('Source changed before copy')
                for block in iter(lambda:f.read(4*1024**2),b''):
                    g.write(block);h.update(block);rec['file_bytes_done']+=len(block);rec['copied_bytes']+=len(block)
                    if time.monotonic()-last>5:save();last=time.monotonic()
                g.flush();os.fsync(g.fileno())
            if identity(src)!=item['identity'] or rec['file_bytes_done']!=item['identity']['bytes']:raise ValueError('Source changed during copy')
            rec.update(stage='hashing_target',file_bytes_done=0);save()
            target_hash=hashlib.sha256();last=time.monotonic()
            with part.open('rb') as f:
                for block in iter(lambda:f.read(4*1024**2),b''):
                    target_hash.update(block);rec['file_bytes_done']+=len(block)
                    if time.monotonic()-last>5:save();last=time.monotonic()
            if target_hash.hexdigest()!=h.hexdigest():raise ValueError('Copied content differs')
            part.rename(dst)
            rec['accepted_files'].append(dict(source=str(src),destination=str(dst),bytes=item['identity']['bytes'],sha256=h.hexdigest(),source_identity=item['identity']))
            save()
        rewrites=[]
        for entry in meta['metadata_files']:
            p=target/entry['path'];raw=p.read_bytes()
            if hashlib.sha256(raw).hexdigest()!=entry['sha256'] or sha(Path(entry['source']))!=entry['sha256']:raise ValueError('Metadata changed')
            new=rebase(raw,maps)
            if str(source).encode() in new:raise ValueError('Unresolved original checkpoint path')
            if new!=raw:
                p.write_bytes(new)
                if sha(p)!=hashlib.sha256(new).hexdigest():raise ValueError('Metadata rewrite differs')
                rewrites.append(dict(path=entry['path'],old_sha256=entry['sha256'],new_sha256=sha(p)))
        if (source/'_lock').exists():raise ValueError('Source acquired runtime lock')
        (target/'.inspection_only.json').rename(target/'.pre_migration_metadata_receipt.json')
        rec.update(stage='migrated_pending_recovery_validation',metadata_rewrites=rewrites,
                   source_metadata_unchanged=True,payload_bytes=plan['expected_bytes'],execution_allowed=False)
        atomic(control/'acceptance.json',rec)
    except BaseException as exc:
        rec.update(stage='failed_preserved',error=repr(exc));raise
    finally:
        save();(control/'owner.lock').unlink()
    return rec

def main():
    p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--control',type=Path,required=True)
    a=p.parse_args();migrate(json.loads(a.plan.read_text()),a.control)
if __name__=='__main__':main()
