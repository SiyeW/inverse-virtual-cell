"""Make a bounded, inspection-only metadata copy. Never copy payloads or run jobs."""
import argparse,hashlib,json,os,re,stat
from pathlib import Path

def clone(source,target):
    source=source.resolve(strict=True);target=target.absolute()
    if target.exists() or target.is_symlink() or target.resolve().is_relative_to(source):
        raise ValueError('Independent unused inspection target required')
    if not (source/'_invocation').is_file() or not (source/'_mrosource').is_file():raise ValueError('Not a recorded pipeline')
    if (source/'_lock').exists():raise ValueError('Source runtime lock remains; inspect ownership first')
    target.mkdir();records=[];aliases=[];count=0;total=0
    def visit(src,dst,ancestors):
        nonlocal count,total
        real=src.resolve(strict=True)
        if not real.is_relative_to(source) or real in ancestors:raise ValueError('Escaping or cyclic directory alias')
        count+=1
        if count>2000:raise ValueError('Directory budget exceeded')
        ancestors=ancestors|{real}
        for item in sorted(src.iterdir()):
            out=dst/item.name
            if item.is_dir():
                if item.name in ['files','tmp','journal','outs']:
                    out.mkdir();continue
                m=re.fullmatch(r'(.+)-u[0-9a-f]{10}',item.name)
                if m:
                    alias=src/m.group(1)
                    if alias.is_symlink() and alias.resolve()==item.resolve():continue
                if item.is_symlink():aliases.append(dict(path=str(out.relative_to(target)),source_target=str(item.resolve().relative_to(source))))
                out.mkdir();visit(item,out,ancestors)
            elif item.name.startswith('_') and item.name not in ['_lock','_pid']:
                if not item.resolve().is_relative_to(source):raise ValueError('Escaping metadata link')
                if not stat.S_ISREG(item.stat().st_mode):raise ValueError('Non-regular metadata refused')
                with item.open('rb') as f:
                    before=os.fstat(f.fileno())
                    if before.st_size>2_000_000:raise ValueError('Metadata file budget exceeded')
                    data=f.read(2_000_001);after=os.fstat(f.fileno())
                if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino) or len(data)!=before.st_size:
                    raise ValueError('Source metadata changed')
                total+=len(data)
                if total>2_000_000:raise ValueError('Total metadata budget exceeded')
                h=hashlib.sha256(data).hexdigest()
                with out.open('xb') as f:f.write(data)
                if hashlib.sha256(out.read_bytes()).hexdigest()!=h:raise ValueError('Metadata copy hash differs')
                records.append(dict(path=str(out.relative_to(target)),source=str(item),bytes=len(data),sha256=h))
    visit(source,target,set())
    if (source/'_lock').exists():raise ValueError('Source acquired runtime lock during inspection copy')
    receipt=dict(source=str(source),target=str(target),execution_allowed=False,payloads_copied=False,
                 output_paths_rewritten=False,metadata_bytes=total,metadata_files=records,aliases_flattened=aliases)
    (target/'.inspection_only.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--target',type=Path,required=True)
    a=p.parse_args();r=clone(a.source,a.target);print(json.dumps({k:v for k,v in r.items() if k not in ['metadata_files','aliases_flattened']}))

if __name__=='__main__':main()
