"""One explicit, guarded resume of a validated migrated pipeline; no automatic retry."""
import argparse,hashlib,json,os,shutil,signal,subprocess,sys,time
from datetime import datetime
from pathlib import Path
from migrate_checkpoint import atomic,sha
import copy_fallback

def check_inputs(items):
    for x in items:
        p=Path(x['path']);s=p.stat();expected=x['identity']
        for field,value in [('bytes',s.st_size),('mtime_ns',s.st_mtime_ns),('device',s.st_dev)]:
            if expected[field]!=value:raise ValueError('Input identity differs: '+str(p)+' '+field)

def validate(config):
    root=Path(config['root']).resolve(strict=True);run=Path(config['run']).resolve(strict=True)
    if not run.is_relative_to(root) or 'data' not in root.parts:raise ValueError('Wrong working root')
    if Path(config['source']).resolve().is_relative_to(root):raise ValueError('Old source overlaps new work')
    for p in [root/'active_cellranger.lock',run/'_lock',root/'control/migration/owner.lock',Path(config['source'])/'_lock']:
        if p.exists():raise ValueError('Existing writer lock: '+str(p))
    acceptance=root/'control/migration/acceptance.json'
    if sha(acceptance)!=config['migration_sha256']:raise ValueError('Migration proof differs')
    migration=json.loads(acceptance.read_text())
    if (migration['stage']!='migrated_pending_recovery_validation' or len(migration['accepted_files'])!=12
        or sum(x['bytes'] for x in migration['accepted_files'])!=14076787910 or migration['target']!=str(run)):
        raise ValueError('Incomplete migration')
    for x in migration['accepted_files']:
        p=Path(x['destination'])
        if not p.is_relative_to(run) or p.stat().st_size!=x['bytes']:raise ValueError('Migrated payload changed')
    meta=json.loads((run/'.pre_migration_metadata_receipt.json').read_text())
    rewrites={x['path']:x['new_sha256'] for x in migration['metadata_rewrites']}
    for x in meta['metadata_files']:
        if sha(run/x['path'])!=rewrites.get(x['path'],x['sha256']):raise ValueError('Checkpoint metadata changed')
    shard=run/config['shard_relative']
    if not (shard/'chnk1/_complete').exists() or (shard/'chnk0/_complete').exists():raise ValueError('Wrong completed/failed siblings')
    if not (shard/'chnk0/_errors').exists():raise ValueError('Expected preserved RNA failure')
    out=json.loads((shard/'chnk1/_outs').read_text())
    if out['total_read_pairs']!=78839829:raise ValueError('Completed CRISPR count differs')
    if any(Path(p).is_symlink() for p in [shard/'chnk0',shard/'chnk1']):raise ValueError('Unexpected aliases')
    if (shard/'chnk0/files').exists() and any((shard/'chnk0/files').iterdir()):raise ValueError('Failed RNA payload should not have migrated')
    cfg=Path(config['csv']);raw=cfg.read_bytes()
    if sha(cfg)!=config['csv_sha256'] or b'force-cells' in raw or b'create-bam,false' not in raw:raise ValueError('Scientific recipe changed')
    check_inputs(config['inputs'])
    view=Path(config['tool']);v=json.loads((view/'compat_view.json').read_text())
    if sha(view/'compat_view.json')!=config['view_sha256']:raise ValueError('Software view proof changed')
    for x in v['files']:
        if x['kind'] in ['verified_copy','fixed_loader_config','local_only_report_copy_policy'] and sha(view/x['path'])!=x['sha256']:
            raise ValueError('Compatibility software changed')
    for p,h in config['engineering_receipts'].items():
        if sha(Path(p))!=h:raise ValueError('Engineering acceptance changed')
    return root,run,shard

def work_bytes(root):
    total=0
    for directory,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=[d for d in dirs if not (Path(directory)/d).is_symlink()]
        for name in files:
            try:total+=(Path(directory)/name).lstat().st_size
            except FileNotFoundError:pass
    return total

def unexpected_completed_reexecution(lines):
    return any('(run:local)' in x and '.MAKE_SHARD.fork0.chnk1.main' in x for x in lines.splitlines())

def run(config):
    root,output,shard=validate(config);policy=config['budgets']
    if any(type(policy.get(k)) is not int or policy[k]<=0 for k in ['localcores','localmem_gb','maximum_work_bytes','minimum_free_bytes']):raise ValueError('Explicit resource limits required')
    available=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024
    if available<(policy['localmem_gb']+16)*1024**3:raise ValueError('Insufficient shared RAM headroom')
    if shutil.disk_usage(root).free<policy['maximum_work_bytes']+policy['minimum_free_bytes']:raise ValueError('Insufficient data space')
    if policy['localcores']>len(os.sched_getaffinity(0)):raise ValueError('CPU budget exceeds affinity')
    control=root/'control/resume';control.mkdir()
    argv=[str(Path(config['tool'])/'cellranger'),'multi','--id','chip02_full_r1','--csv',config['csv'],
          '--output-dir',str(output),'--localcores',str(policy['localcores']),'--localmem',str(policy['localmem_gb']),'--disable-ui']
    env=dict(os.environ,MRO_UNIQUIFIED_DIRECTORIES='disable',MRO_CIFS_COMPAT='1',MRO_CIFS_SAME_MODE_ROOT=str(root))
    env.pop('MRO_FULLSTAGERESET',None)
    cmd=[sys.executable,str(Path(copy_fallback.__file__)),'--']+argv
    rec=dict(stage='launching',pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19],
             command=cmd,config=config,started_at=datetime.now().astimezone().isoformat(),reports_complete=False,
             old_source_preserved=True,failed_rna_chunk_reexecution_expected=True,completed_crispr_chunk_reexecution_allowed=False)
    lock=root/'active_cellranger.lock'
    with lock.open('x') as f:json.dump(dict(pid=rec['pid'],uid=rec['uid'],startticks=rec['startticks']),f)
    atomic(control/'launch.json',rec)
    proc=None;started=time.monotonic();last=0;read_offset=0
    def save():
        rec.update(updated_at=datetime.now().astimezone().isoformat(),elapsed_seconds=time.monotonic()-started);atomic(control/'status.json',rec)
    try:
        with (control/'stdout').open('x') as out,(control/'stderr').open('x') as err:
            proc=subprocess.Popen(cmd,env=env,cwd=root,stdout=out,stderr=err,start_new_session=True)
            rec.update(stage='running',pipeline_pid=proc.pid);save()
            while proc.poll() is None:
                with (control/'stdout').open() as log:
                    log.seek(read_offset);new=log.read(256000);read_offset=log.tell()
                if unexpected_completed_reexecution(new):raise RuntimeError('Unexpected completed CRISPR chunk reexecution')
                if time.monotonic()-last>=60 or not rec.get('work_bytes'):
                    rec['work_bytes']=work_bytes(output);last=time.monotonic()
                rec['free_bytes']=shutil.disk_usage(root).free
                if rec['work_bytes']>policy['maximum_work_bytes'] or rec['free_bytes']<policy['minimum_free_bytes']:
                    raise RuntimeError('Disk guard reached; preserve and stop')
                save();time.sleep(10)
        rec['exit_code']=proc.returncode
        if proc.returncode:raise RuntimeError('Pipeline exited nonzero; no automatic retry')
        rec['stage']='cellranger_complete_pending_acceptance'
    except BaseException as exc:
        if proc is not None and proc.poll() is None:
            os.killpg(proc.pid,signal.SIGINT)
            try:proc.wait(timeout=60)
            except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=30)
        rec.update(stage='stopped_preserved',error=repr(exc));raise
    finally:save();lock.unlink()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()))
