"""Create an isolated software view; never overwrite an installed component."""
import argparse, hashlib, json, os, shutil
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def stage(source, runtime, target):
    source=source.resolve(strict=True);runtime=runtime.resolve(strict=True)
    target=target.absolute()
    if target.exists() or target.is_symlink() or target.is_relative_to(source):
        raise ValueError('An independent unused target is required')
    for rel in ['bin/cellranger','external/martian/bin/mrjob']:
        if not (source/rel).is_file():raise ValueError('Incomplete original installation')
    if not runtime.is_file():raise ValueError('Missing tested runtime')
    target.mkdir()
    custom={'bin/cellranger':source/'bin/cellranger','external/martian/bin/mrp':runtime}
    expand={'bin','external','external/martian','external/martian/bin'}
    receipts=[]
    def visit(src,dst,prefix=''):
        for child in sorted(src.iterdir()):
            rel=(prefix+'/'+child.name).strip('/');out=dst/child.name
            if rel in custom:
                origin=custom[rel];h=digest(origin);shutil.copyfile(origin,out)
                if digest(out)!=h:raise ValueError('Executable copy hash mismatch')
                out.chmod(0o755)
                receipts.append(dict(path=rel,source=str(origin),sha256=h,kind='verified_copy'))
            elif rel in expand:
                out.mkdir();visit(child,out,rel)
            else:
                os.symlink(str(child),out,target_is_directory=child.is_dir())
                receipts.append(dict(path=rel,source=str(child),kind='source_reference_not_copied'))
    visit(source,target)
    # The launcher resolves its executable path; point this entry to our copied
    # launcher rather than following the original installation's entry symlink.
    entry=target/'cellranger'
    if not entry.is_symlink():raise ValueError('Unexpected launcher layout')
    entry.unlink();entry.symlink_to('bin/cellranger')
    record=dict(source=str(source),target=str(target),original_modified=False,
                scientific_binaries_modified=False,files=receipts,experimental=True)
    (target/'compat_view.json').write_text(json.dumps(record,indent=2)+'\n')
    return record

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--runtime',type=Path,required=True);p.add_argument('--target',type=Path,required=True)
    a=p.parse_args();print(json.dumps(stage(a.source,a.runtime,a.target)))

if __name__=='__main__':main()
