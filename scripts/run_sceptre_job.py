"""Own one official SCEPTRE calibration process with durable bounded status."""
from pathlib import Path
import argparse,hashlib,json,os,socket,subprocess,time

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for part in iter(lambda:f.read(1048576),b''):h.update(part)
    return h.hexdigest()

def startticks(pid):
    s=Path(f'/proc/{pid}/stat').read_text()
    return int(s[s.rfind(')')+2:].split()[19])

def run(cfg):
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
    native=Path(cfg['native']);output=Path(cfg['output'])
    assert native.is_dir() and not (native/'run.exit').exists()
    lock=native/'active.lock';fd=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL)
    os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=startticks(os.getpid()),
        job='official_sceptre_'+cfg.get('mode','NTC_calibration'),commit=cfg['commit'])).encode());os.close(fd)
    child=None
    def state(stage,**kw):
        p=native/'status.json';tmp=p.with_suffix('.receiving')
        tmp.write_text(json.dumps(dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
            method=cfg['mechanism'],scientific_parent_complete=False,**kw),indent=2));tmp.replace(p)
    try:
        state('verify_accepted_prepared_counts_and_published_sources')
        for path,want in cfg['sha256'].items():
            if sha(path)!=want:raise ValueError('Frozen source differs: '+path)
        proof=json.loads(Path(cfg['input_proof']).read_text())
        expected=cfg.get('expected',dict(records=71664,NTC_records=19342,NTC_guides=106,candidate_guides=678))
        assert proof['records']==expected['records'] and proof['response_RNAs']==2000 and proof['NTC_records']==expected['NTC_records']
        assert proof['recovery']['roundtrip_bit_exact'] and proof['recovery']['adjacent_integer_uniqueness']
        for name,item in proof['files'].items():
            if sha(Path(cfg['input'])/name)!=item['sha256']:raise ValueError('Prepared count input differs: '+name)
        assert not output.exists();output.mkdir()
        env=os.environ|json.loads(Path(cfg['runtime_environment']).read_text())|{
            'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2','MKL_NUM_THREADS':'2','R_MAX_NUM_DLLS':'200'}
        screen=cfg.get('mode')=='guide_screen'
        command=(['/usr/local/bin/Rscript','--vanilla',cfg['R_script'],cfg['summary_R_script'],cfg['input'],cfg['calibration_object'],str(output),str(cfg['seed'])] if screen else
            ['/usr/local/bin/Rscript','--vanilla',cfg['R_script'],'--run',cfg['input'],str(output),
             cfg['mechanism'],str(cfg['pairs']),str(cfg['seed']),
             ','.join(str(expected[k]) for k in ['records','NTC_records','NTC_guides','candidate_guides'])])
        begin=time.monotonic()
        with (native/'R.stdout').open('w') as out,(native/'R.stderr').open('w') as err:
            child=subprocess.Popen(command,env=env,stdout=out,stderr=err)
            identity=dict(pid=child.pid,uid=os.getuid(),startticks=startticks(child.pid))
            while child.poll() is None:
                phase=output/'stage.txt'
                state('running_official_sceptre',phase=phase.read_text().strip() if phase.exists() else 'R_import',
                    R_owner=identity,elapsed_seconds=round(time.monotonic()-begin,1))
                if time.monotonic()-begin>cfg['science_guard_seconds']:
                    child.terminate()
                    try:child.wait(timeout=30)
                    except subprocess.TimeoutExpired:child.kill();child.wait(timeout=30)
                    raise TimeoutError('Scientific protection reached; completed input cache is retained')
                time.sleep(15)
        if child.returncode:
            state('failed_pending_recovery',exit=child.returncode)
            return child.returncode
        required=['calibration_summary.csv','NTC_expression_tests.csv','cohort_contract.tsv','sessionInfo.txt']
        targeting=False
        if screen:
            required=['NTC_guide_summary.csv','gate_status.tsv','sessionInfo.txt']
            targeting=(output/'guide_screen_summary.csv').exists()
            if targeting:required+=['guide_screen_summary.csv','same_target_pairs.csv','screen_summary.csv']
            if targeting and cfg.get('energy_config'):
                state('compute_fixed_PC_effect_sizes_after_screen')
                with (native/'energy.stdout').open('w') as out,(native/'energy.stderr').open('w') as err:
                    energy=subprocess.run([cfg['scientific_python'],cfg['energy_script'],'--config',cfg['energy_config']],
                        env=os.environ|{'OPENBLAS_NUM_THREADS':'2','OMP_NUM_THREADS':'2','NUMBA_NUM_THREADS':'2','NUMEXPR_NUM_THREADS':'2'},
                        stdout=out,stderr=err,timeout=600)
                if energy.returncode:raise RuntimeError('Effect-size export failed; completed expression tests are retained')
                required+=['energy_by_guide.csv','energy_proof.json']
        for name in required:assert (output/name).is_file()
        receipt=dict(commit=cfg['commit'],official_sceptre_tag='v0.10.3',
            official_sceptre_commit='81d32c6e6e94a9aca84bccbe484ec869ef769087',
            source_input_proof_sha256=sha(cfg['input_proof']),no_targeting_tests_performed=not targeting,
            input_scope=proof['input_scope'],full_frozen_NTC_reference=True,
            additional_cell_sampling_or_QC=False,
            files={name:dict(bytes=(output/name).stat().st_size,sha256=sha(output/name)) for name in required})
        (output/'result_receipt.json').write_text(json.dumps(receipt,indent=2))
        state('complete_pending_scientific_review',elapsed_seconds=round(time.monotonic()-begin,1))
        return 0
    except BaseException as exc:
        if child is not None and child.poll() is None:
            child.terminate()
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:child.kill();child.wait(timeout=30)
        state('failed_pending_recovery',error=type(exc).__name__+': '+str(exc))
        raise
    finally:lock.unlink(missing_ok=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args()
    raise SystemExit(run(json.loads(a.config.read_text())))
