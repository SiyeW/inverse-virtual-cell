"""Exercise a completed and a failed sibling chunk across a new filesystem root."""
import argparse,json,os,subprocess,time
from pathlib import Path
from run_resume_fixture import stop
from migrate_checkpoint import migrate,identity,sha

def run(old,new,runtime,python_bin):
    if old.exists() or new.exists():raise ValueError('Unused fixture roots required')
    old.mkdir(parents=True);new.mkdir(parents=True)
    fixture=Path(__file__).with_name('relocation_fixture').resolve()
    marker=old/'once';release=new/'release'
    call=old/'call.mro';call.write_text('@include "pipeline.mro"\ncall CIFS_MULTI_RECOVERY(marker='+json.dumps(str(marker))+', release='+json.dumps(str(release))+',)\n')
    env=dict(os.environ,PATH=str(python_bin)+':'+os.environ['PATH'],MROPATH=str(fixture))
    env.pop('MRO_UNIQUIFIED_DIRECTORIES',None);env.pop('MRO_CIFS_COMPAT',None)
    source=old/'run';target=new/'run'
    def cmd(root,inv):return [str(runtime),str(inv),'recovery_test','--psdir='+str(root),'--disable-ui','--localcores=2','--localmem=3']
    first=second=None;rec=dict(complete=False,science=False)
    try:
        with (old/'stdout').open('w') as out,(old/'stderr').open('w') as err:
            first=subprocess.Popen(cmd(source,call),env=env,cwd=old,stdout=out,stderr=err,start_new_session=True)
            chunk=source/'CIFS_MULTI_RECOVERY/SPLIT_WORK/fork0/chnk1'
            deadline=time.monotonic()+90
            while not (chunk/'_complete').exists():
                if first.poll() is not None or time.monotonic()>deadline:raise RuntimeError('Fixture did not complete sibling')
                time.sleep(.5)
            stop(first)
        for _ in range(40):
            try:os.killpg(first.pid,0)
            except ProcessLookupError:break
            time.sleep(.25)
        else:raise RuntimeError('Old fixture group remains')
        outs=json.loads((chunk/'_outs').read_text());payload=Path(outs['piece']).resolve()
        plan=dict(source=str(source),target=str(target),expected_files=1,expected_bytes=payload.stat().st_size,
                  files=[dict(path=str(payload),identity=identity(payload))])
        r=migrate(plan,new/'migration')
        old_hashes={str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file() and not p.is_symlink()}
        release.write_text('release failed sibling\n')
        env.update(MRO_UNIQUIFIED_DIRECTORIES='disable',MRO_CIFS_COMPAT='1')
        with (new/'stdout').open('w') as out,(new/'stderr').open('w') as err:
            second=subprocess.Popen(cmd(target,target/'_invocation'),env=env,cwd=new,stdout=out,stderr=err,start_new_session=True)
            second.wait(timeout=120)
        if second.returncode:raise RuntimeError('Relocated fixture failed')
        results=list((target/'outs').iterdir());assert len(results)==1 and results[0].read_bytes()==b'41\n'
        assert all(sha(source/p)==h for p,h in old_hashes.items())
        logs=(new/'stdout').read_text();assert '(reset-partial)' in logs
        assert not any('(run:local)' in line and 'chnk1' in line for line in logs.splitlines())
        rec.update(complete=True,completed_sibling_not_rerun=True,source_unchanged=True,final_output='41',migration=r)
    finally:
        for p in [first,second]:
            if p is not None:stop(p)
        (new/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    return rec

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',type=Path,required=True);p.add_argument('--new',type=Path,required=True);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--python-bin',type=Path,required=True)
    a=p.parse_args();print(json.dumps(run(a.old,a.new,a.runtime,a.python_bin)))
