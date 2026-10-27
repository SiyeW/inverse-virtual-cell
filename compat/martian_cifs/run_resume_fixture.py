"""Bounded synthetic interruption test; cannot accept a pre-existing run."""
import argparse,hashlib,json,os,signal,subprocess,time
from itertools import islice
from pathlib import Path

def stop(proc):
    if proc.poll() is None:
        os.killpg(proc.pid,signal.SIGINT)
        try:proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=10)

def proc_identity(pid):
    p=Path('/proc')/str(pid)
    try:
        if p.stat().st_uid!=os.getuid():return None
        with (p/'stat').open() as f:fields=f.read(8192).rsplit(')',1)[1].split()
        return dict(pid=pid,uid=os.getuid(),startticks=fields[19],state=fields[0],pgrp=int(fields[2]))
    except FileNotFoundError:return None

def owned_tree(pid):
    queue=[pid];rows=[];seen=set()
    while queue:
        pid=queue.pop(0)
        if pid in seen:continue
        if len(seen)>=32:raise RuntimeError('Synthetic process tree exceeds explicit bound')
        seen.add(pid);row=proc_identity(pid)
        if not row:continue
        rows.append(row)
        try:
            tasks=list(islice((Path('/proc')/str(pid)/'task').iterdir(),257))
            if len(tasks)>256:raise RuntimeError('Synthetic thread bound exceeded')
            for task in tasks:
                try:
                    with (task/'children').open() as stream:raw=stream.read(8193)
                    if len(raw)>8192:raise RuntimeError('Synthetic child list exceeds bound')
                    queue.extend(int(x) for x in raw.split())
                except FileNotFoundError:pass
        except FileNotFoundError:pass
    return rows

def still_live(known,reader=proc_identity):
    live=[]
    for old in known:
        row=reader(old['pid'])
        if row and (row['uid'],row['startticks'])==(old['uid'],old['startticks']) and row['state']!='Z':
            live.append(row)
    return live

def run(case,runtime,python_bin,phase_timeout=90):
    if not 30<=phase_timeout<=300:raise ValueError('Bounded fixture phase timeout required')
    case=case.absolute();runtime=runtime.resolve(strict=True);python_bin=python_bin.resolve(strict=True)
    if case.exists() or case.is_symlink() or 'provenance' not in case.parts:
        raise ValueError('An unused project provenance test directory is required')
    case.mkdir();fixture=Path(__file__).with_name('resume_fixture').resolve()
    marker=case/'make_once.marker';started=case/'wait_started.marker'
    invocation=case/'call.mro'
    invocation.write_text('@include "pipeline.mro"\ncall CIFS_RESUME_TEST(marker='+json.dumps(str(marker))+', started='+json.dumps(str(started))+',)\n')
    env=dict(os.environ,PATH=str(python_bin)+':'+os.environ['PATH'],MROPATH=str(fixture),
             MRO_UNIQUIFIED_DIRECTORIES='disable',MRO_CIFS_COMPAT='1')
    cmd=[str(runtime),str(invocation),'resume_fixture','--psdir='+str(case/'run'),
         '--disable-ui','--localcores=1','--localmem=2']
    rec=dict(command=cmd,scientific_result=False,laboratory_data_read=False)
    first=None;second=None
    try:
        with (case/'first.stdout').open('wb') as out,(case/'first.stderr').open('wb') as err:
            first=subprocess.Popen(cmd,cwd=case,env=env,stdout=out,stderr=err,start_new_session=True)
            deadline=time.monotonic()+phase_timeout
            while not started.exists():
                if first.poll() is not None:raise RuntimeError('Initial synthetic run ended before planned interruption')
                if time.monotonic()>deadline:raise TimeoutError('No waiting-stage marker')
                time.sleep(.5)
            assert marker.read_bytes()==b'executed once\n'
            rec['completed_step_marker_sha256']=hashlib.sha256(marker.read_bytes()).hexdigest()
            known=owned_tree(first.pid)
            rec['first_pid']=first.pid;rec['known_owned_processes']=known
            stop(first);rec['first_exit_code']=first.returncode
            if first.returncode==0:raise RuntimeError('Run completed before interruption')
        # No overlapping writer: both the driver and its process group must end.
        for _ in range(80):
            if not still_live(known):break
            time.sleep(.25)
        else:raise RuntimeError('Prior positively identified synthetic writer remains live')
        rec['old_writers_exited']=True
        # Zombies cannot write; killpg(..., 0) can nevertheless report their
        # group as present until PID1 reaps them. Do not kill unrelated services.
        try:os.killpg(first.pid,0)
        except ProcessLookupError:rec['old_process_group_exited']=True
        else:rec['old_process_group_exited']=False
        with (case/'resume.stdout').open('wb') as out,(case/'resume.stderr').open('wb') as err:
            second=subprocess.Popen(cmd,cwd=case,env=env,stdout=out,stderr=err,start_new_session=True)
            rec['resume_exit_code']=second.wait(timeout=phase_timeout)
        if second.returncode:raise RuntimeError('Synthetic resume failed')
        assert marker.read_bytes()==b'executed once\n'
        results=[p for p in (case/'run/outs').iterdir() if p.is_file()]
        assert len(results)==1 and results[0].read_bytes()==b'42\n'
        rec.update(complete=True,completed_step_not_reexecuted=True,
                   output_sha256=hashlib.sha256(results[0].read_bytes()).hexdigest())
    except Exception as exc:
        rec.update(complete=False,error=repr(exc));raise
    finally:
        for proc in [first,second]:
            if proc is not None:stop(proc)
        (case/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    return rec

def main():
    p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True)
    p.add_argument('--runtime',type=Path,required=True);p.add_argument('--python-bin',type=Path,required=True)
    p.add_argument('--phase-timeout',type=int,default=90)
    a=p.parse_args();print(json.dumps(run(a.case,a.runtime,a.python_bin,a.phase_timeout)))

if __name__=='__main__':main()
