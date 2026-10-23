"""Bounded synthetic interruption test; cannot accept a pre-existing run."""
import argparse,hashlib,json,os,signal,subprocess,time
from pathlib import Path

def stop(proc):
    if proc.poll() is None:
        os.killpg(proc.pid,signal.SIGINT)
        try:proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=10)

def run(case,runtime,python_bin):
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
            deadline=time.monotonic()+90
            while not started.exists():
                if first.poll() is not None:raise RuntimeError('Initial synthetic run ended before planned interruption')
                if time.monotonic()>deadline:raise TimeoutError('No waiting-stage marker')
                time.sleep(.5)
            assert marker.read_bytes()==b'executed once\n'
            rec['completed_step_marker_sha256']=hashlib.sha256(marker.read_bytes()).hexdigest()
            stop(first);rec['first_exit_code']=first.returncode
            if first.returncode==0:raise RuntimeError('Run completed before interruption')
        # No overlapping writer: both the driver and its process group must end.
        for _ in range(20):
            try:os.killpg(first.pid,0)
            except ProcessLookupError:break
            time.sleep(.25)
        else:raise RuntimeError('Prior owned process group still exists')
        rec['old_process_group_exited']=True
        with (case/'resume.stdout').open('wb') as out,(case/'resume.stderr').open('wb') as err:
            second=subprocess.Popen(cmd,cwd=case,env=env,stdout=out,stderr=err,start_new_session=True)
            rec['resume_exit_code']=second.wait(timeout=90)
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
    a=p.parse_args();print(json.dumps(run(a.case,a.runtime,a.python_bin)))

if __name__=='__main__':main()
