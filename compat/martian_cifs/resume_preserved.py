"""One explicit resume after a frozen metadata inventory; never retry blindly."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime

from migrate_checkpoint import atomic, sha
from resume_migrated import check_inputs, work_bytes
import copy_fallback


def owned_relative(run, value):
    relative=Path(value)
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Nonrelative preserved path')
    return run/relative


def preserved_metadata_matches(run, inventory):
    """Read only frozen small metadata and stat declared completed outputs."""
    run = Path(run)
    for node in inventory['preserved_complete_nodes']:
        if not node['complete']:
            raise ValueError('Noncompleted node in preserved set')
        for item in node['metadata'].values():
            p = owned_relative(run,item['relative'])
            if p.is_symlink() or p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
                raise ValueError('Completed metadata changed: ' + str(p))
        for item in node['outputs']:
            p = owned_relative(run,item['relative'])
            if not p.is_relative_to(run) or p.is_symlink():
                raise ValueError('Output path outside preserved run')
            current = p.stat()
            if (current.st_size, current.st_mtime_ns, current.st_dev) != (item['bytes'], item['mtime_ns'], item['device']):
                raise ValueError('Completed output identity changed: ' + str(p))
    for key in ['UUID', 'invocation']:
        item = inventory[key]
        if sha(owned_relative(run,item['relative'])) != item['sha256']:
            raise ValueError('Pipestance identity changed: ' + key)
    for item in inventory.get('required_resume_inputs', []):
        p=owned_relative(run,item['relative']);current=p.stat()
        if (current.st_size,current.st_mtime_ns,current.st_dev)!=(item['bytes'],item['mtime_ns'],item['device']):
            raise ValueError('Required remaining-work input changed: '+str(p))


def completed_node_reexecution(text, run_name, inventory):
    forbidden = ['.' + node['relative'].replace('/', '.') + '.'
                 for node in inventory['preserved_complete_nodes']]
    for line in text.splitlines():
        if '(run:local)' not in line:
            continue
        if any(key in line for key in forbidden):
            return line
        if '.MAKE_SHARD.' in line or '.BARCODE_CORRECTION.' in line:
            return line
        if any(fragment in line for fragment in inventory.get('forbidden_stage_fragments', [])):
            return line
    return None


def validate(config):
    root=Path(config['root']).resolve(strict=True)
    run=Path(config['run']).resolve(strict=True)
    if not run.is_relative_to(root) or 'data' not in root.parts or run.name != 'chip02_full_r1':
        raise ValueError('Unexpected explicit production root')
    if config.get('user_authorization') != 'repair_then_resume_preserving_completed_chip02':
        raise ValueError('Explicit checkpoint-preserving recovery required')
    for p in [root/'active_cellranger.lock',run/'_lock']:
        if p.exists():raise ValueError('Existing output owner: '+str(p))
    audit=Path(config['inventory'])
    if sha(audit)!=config['inventory_sha256']:raise ValueError('Recovery inventory identity differs')
    inventory=json.loads(audit.read_text())
    if inventory['run']!=str(run) or inventory['source_payload_read']:
        raise ValueError('Wrong recovery inventory')
    if sorted(inventory['complete_chunks'])!=config['preserve_alignment_chunks']:
        raise ValueError('Completed chunk set changed')
    if sorted(inventory['incomplete_chunks'])!=config['redo_alignment_chunks']:
        raise ValueError('Incomplete chunk set changed')
    preserved_metadata_matches(run,inventory)
    check_inputs(config['inputs'])
    csv=Path(config['csv']);raw=csv.read_bytes()
    if sha(csv)!=config['csv_sha256'] or b'force-cells' in raw or b'create-bam,false' not in raw:
        raise ValueError('Scientific recipe changed')
    view=Path(config['tool'])
    if sha(view/'compat_view.json')!=config['view_sha256']:
        raise ValueError('Private software view identity changed')
    software=json.loads((view/'compat_view.json').read_text())
    for item in software['files']:
        if item['kind'] in ['verified_copy','fixed_loader_config','local_only_report_copy_policy','local_only_directory_copy_policy']:
            if sha(view/item['path'])!=item['sha256']:raise ValueError('Software copy changed')
    if sha(Path(config['repair_receipt']))!=config['repair_receipt_sha256']:
        raise ValueError('Repair tests receipt changed')
    control=root/config['control_relative']
    if not control.is_relative_to(root/'control') or control.exists():
        raise ValueError('A fresh control directory is required')
    policy=config['budgets']
    if policy!={'localcores':8,'localmem_gb':64,'maximum_work_bytes':1000000000000,'minimum_free_bytes':100000000000}:
        raise ValueError('Existing production resource limits must remain unchanged')
    return root,run,control,inventory


def run(config):
    root,output,control,inventory=validate(config)
    policy=config['budgets']
    available=int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024
    if available<(policy['localmem_gb']+16)*1024**3:raise ValueError('Insufficient shared RAM headroom')
    if shutil.disk_usage(root).free < policy['maximum_work_bytes']+policy['minimum_free_bytes']:
        raise ValueError('Insufficient existing disk safety headroom')
    if policy['localcores']>len(os.sched_getaffinity(0)):raise ValueError('CPU affinity smaller than existing budget')
    control.mkdir()
    cmd=[sys.executable,str(Path(copy_fallback.__file__)),'--',str(Path(config['tool'])/'cellranger'),
         'multi','--id','chip02_full_r1','--csv',config['csv'],'--output-dir',str(output),
         '--localcores','8','--localmem','64','--disable-ui']
    env=dict(os.environ,MRO_UNIQUIFIED_DIRECTORIES='disable',MRO_CIFS_COMPAT='1',MRO_CIFS_SAME_MODE_ROOT=str(root))
    env.pop('MRO_FULLSTAGERESET',None)
    record=dict(stage='launching',pid=os.getpid(),uid=os.getuid(),
        startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19],command=cmd,
        config=config,started_at=datetime.now().astimezone().isoformat(),reports_complete=False,
        old_source_preserved=True,completed_alignment_chunks=inventory['complete_chunks'],
        incomplete_alignment_chunks=inventory['incomplete_chunks'])
    lock=root/'active_cellranger.lock'
    with lock.open('x') as stream:json.dump({k:record[k] for k in ['pid','uid','startticks']},stream)
    atomic(control/'launch.json',record)
    process=None;started=time.monotonic();last_measure=0;offset=0
    def save():
        record.update(updated_at=datetime.now().astimezone().isoformat(),elapsed_seconds=time.monotonic()-started)
        atomic(control/'status.json',record)
    try:
        with (control/'stdout').open('x') as out,(control/'stderr').open('x') as err:
            process=subprocess.Popen(cmd,env=env,cwd=root,stdout=out,stderr=err,start_new_session=True)
            record.update(stage='running',pipeline_pid=process.pid);save()
            while process.poll() is None:
                with (control/'stdout').open() as stream:
                    stream.seek(offset);text=stream.read(256000);offset=stream.tell()
                event=completed_node_reexecution(text,output.name,inventory)
                if event:raise RuntimeError('Completed-node reexecution refused: '+event)
                if time.monotonic()-last_measure>=60 or not record.get('work_bytes'):
                    record['work_bytes']=work_bytes(output);last_measure=time.monotonic()
                record['free_bytes']=shutil.disk_usage(root).free
                if record['work_bytes']>policy['maximum_work_bytes'] or record['free_bytes']<policy['minimum_free_bytes']:
                    raise RuntimeError('Existing disk safety limit reached; preserve and stop')
                save();time.sleep(10)
        record['exit_code']=process.returncode
        if process.returncode:raise RuntimeError('Resumed pipeline exited nonzero; evidence preserved')
        record['stage']='cellranger_complete_pending_acceptance'
    except BaseException as exc:
        if process is not None and process.poll() is None:
            os.killpg(process.pid,signal.SIGINT)
            try:process.wait(timeout=60)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGTERM);process.wait(timeout=30)
        record.update(stage='stopped_preserved',error=repr(exc));raise
    finally:
        save();lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True)
    args=parser.parse_args();run(json.loads(args.config.read_text()))
