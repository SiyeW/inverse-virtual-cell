"""Explicit one-chip Cell Ranger execution; no automatic retry or next-chip launch."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

import manage_core_study as study


def limits(runtime):
    keys = ['localcores', 'localmem_gb', 'reserved_workspace_bytes', 'minimum_free_bytes']
    if any(type(runtime.get(k)) is not int or runtime[k] <= 0 for k in keys):
        raise ValueError('Explicit positive CPU, RAM and disk budgets required')
    return {k:runtime[k] for k in keys}


def owned_bytes(root):
    total = 0
    if root.exists():
        for directory, dirs, files in os.walk(root, followlinks=False):
            dirs[:] = [d for d in dirs if not (Path(directory)/d).is_symlink()]
            for name in files:
                try:total += (Path(directory)/name).lstat().st_size
                except FileNotFoundError:pass  # a completed stage may clean its own transient file
    return total


def budget_ok(free, used, policy):
    return free >= policy['minimum_free_bytes'] and used <= policy['reserved_workspace_bytes']


def available_memory():
    return int(next(x for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')).split()[1])*1024


def process_identity():
    return dict(pid=os.getpid(), uid=os.getuid(), startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19])


def prepare_plan(config, chip_id, attempt):
    study.named(chip_id);study.named(attempt)
    chips = [c for c in config['chips'] if c['chip_id'] == chip_id]
    if len(chips) != 1:
        raise ValueError('Select one recorded chip')
    runtime = config['runtime'];policy = limits(runtime)
    root = study.absolute(runtime['work_root'])
    run_id = chip_id + '_' + attempt
    if len(run_id) > 64:
        raise ValueError('Cell Ranger run ID exceeds 64 characters')
    output, control = root/'outputs'/run_id, root/'control'/run_id
    for path in [output, control]:
        if not path.resolve().is_relative_to(root) or path.exists():
            raise ValueError('Output/control must be an unused owned path; no automatic resume')
    csv_path = study.absolute(config['study_root'])/'configs/cellranger'/(chip_id+'.csv')
    expected = study.cellranger_csv(config, chips[0]).encode('utf8')
    if csv_path.read_bytes() != expected:
        raise ValueError('Rendered CSV differs from frozen study settings')
    if b'force-cells' in expected or b'create-bam,false\n' not in expected:
        raise ValueError('Formal runs cannot inherit forced engineering cells or enable BAM')
    argv = [runtime['cellranger'], 'multi', '--id', run_id, '--csv', str(csv_path),
            '--localcores', str(policy['localcores']), '--localmem', str(policy['localmem_gb']),
            '--disable-ui', '--output-dir', str(output)]
    return dict(root=root, output=output, control=control, argv=argv, policy=policy,
                csv_sha256=study.digest(csv_path), run_id=run_id)


def check_authority(config):
    receipt = study.checked_metadata(config['guide_counting_receipt'])
    authority = study.file_spec(config['authoritative_feature_reference'])
    if study.digest(authority) != config['authoritative_feature_reference']['sha256']:
        raise ValueError('Original guide/target authority differs')
    if (receipt.get('capture_fields_unchanged') is not True or receipt.get('target_and_control_relabeling') is not False
        or receipt.get('source_sha256') != config['authoritative_feature_reference']['sha256']
        or receipt['files']['guide_counting_features.csv']['sha256'] != config['feature_reference']['sha256']):
        raise ValueError('Counting reference is not bound to unchanged original guide definitions')
    prefix = study.checked_metadata(config['runtime']['engineering_status'])
    if prefix.get('stage') != 'real_library_prefix_engineering_passed' or prefix.get('scientific_result') is not False:
        raise ValueError('Required bounded real-library engineering check has not passed')
    argv = prefix.get('argv', [])
    if not argv or argv[0] != config['runtime']['cellranger']:
        raise ValueError('Engineering check used a different Cell Ranger executable')
    for flag, value in [('--transcriptome', config['reference']['path']), ('--feature-ref', config['feature_reference']['path'])]:
        if flag not in argv or argv.index(flag)+1 >= len(argv) or argv[argv.index(flag)+1] != value:
            raise ValueError('Engineering check used a different reference recipe')
    probe = study.checked_metadata(config['runtime']['workspace_probe'])
    root = study.absolute(config['runtime']['work_root'])
    if (probe.get('complete') is not True or probe.get('work_root') != str(root)
        or probe.get('owner_uid') != os.getuid() or root.stat().st_uid != os.getuid()
        or probe.get('device') != root.stat().st_dev):
        raise ValueError('Working filesystem differs from the owned accepted probe')
    return dict(guide_counting=receipt, engineering=prefix, workspace=probe)


def run(study_dir, chip_id, attempt):
    config = study.load_initialized(study_dir)
    commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=study.ROOT, text=True).strip()
    if commit != config['code_commit'] or subprocess.check_output(['git','status','--porcelain'], cwd=study.ROOT, text=True).strip():
        raise ValueError('Use the fixed clean pushed execution revision')
    readiness = study.check(config)
    if readiness['blockers']:
        raise ValueError('Study is not ready: '+', '.join(readiness['blockers']))
    check_authority(config)
    plan = prepare_plan(config, chip_id, attempt)
    policy, root = plan['policy'], plan['root']
    free = shutil.disk_usage(root).free
    if free < policy['reserved_workspace_bytes'] + policy['minimum_free_bytes']:
        raise ValueError('Free space cannot cover the declared work budget plus safety floor')
    available = available_memory()
    if available < (policy['localmem_gb'] + 16) * 1024**3:
        raise ValueError('Insufficient currently available RAM with safety headroom')
    if policy['localcores'] > len(os.sched_getaffinity(0)):
        raise ValueError('CPU budget exceeds available affinity')
    lock = root/'active_cellranger.lock'
    identity = process_identity()
    with lock.open('x') as f:json.dump(identity,f)
    control = plan['control'];proc = None;created_control = False
    rec = dict(identity, stage='launching', chip_id=chip_id, attempt=attempt, code_commit=commit,
               argv=plan['argv'], output=str(plan['output']), study=str(study_dir),
               started_at=datetime.now().astimezone().isoformat(), csv_sha256=plan['csv_sha256'],
               budgets=policy, quota_is_exclusively_reserved=False, peak_work_bytes=0,
               phase='cellranger_only_not_core_report', all_reports_complete=False)
    start = time.monotonic();last_size = 0.
    def save():
        rec['updated_at'] = datetime.now().astimezone().isoformat();rec['elapsed_seconds'] = time.monotonic()-start
        temp = control/'status.tmp';temp.write_text(json.dumps(rec,indent=2)+'\n');temp.replace(control/'status.json')
    try:
        control.mkdir()
        created_control = True
        (control/'launch_input_inventory.json').write_text(json.dumps(readiness['inventory'],indent=2)+'\n')
        with (control/'cellranger.stdout').open('w') as out, (control/'cellranger.stderr').open('w') as err:
            proc = subprocess.Popen(plan['argv'], cwd=root, stdout=out, stderr=err, start_new_session=True)
            rec.update(stage='running', pipeline_pid=proc.pid);save()
            while proc.poll() is None:
                free = shutil.disk_usage(root).free
                if time.monotonic()-last_size >= 60:
                    rec['work_bytes'] = owned_bytes(plan['output']);last_size = time.monotonic()
                    rec['peak_work_bytes'] = max(rec['peak_work_bytes'],rec['work_bytes'])
                rec['free_bytes'] = free
                if not budget_ok(free,rec.get('work_bytes',0),policy):
                    raise RuntimeError('Disk guard reached; interrupt owned job and preserve checkpoint')
                save();time.sleep(10)
        rec['exit_code'] = proc.returncode
        if proc.returncode:
            raise RuntimeError('Cell Ranger exited nonzero; preserve output, no automatic retry')
        if list((plan['output']/'outs').rglob('*.bam')):
            raise RuntimeError('Unexpected BAM output; hold acceptance')
        rec['stage'] = 'cellranger_complete_pending_acceptance'
    except Exception as exc:
        if proc is not None and proc.poll() is None:
            try:os.killpg(proc.pid,signal.SIGINT)
            except ProcessLookupError:pass
            try:proc.wait(timeout=60)
            except subprocess.TimeoutExpired:
                try:os.killpg(proc.pid,signal.SIGTERM)
                except ProcessLookupError:pass
                proc.wait(timeout=30)
        rec.update(stage='stopped_preserved',error=repr(exc))
        raise
    finally:
        if created_control:save()
        lock.unlink()
    return rec


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study',required=True,type=Path);p.add_argument('--chip',required=True);p.add_argument('--attempt',required=True)
    a = p.parse_args();print(json.dumps(run(a.study,a.chip,a.attempt)))


if __name__ == '__main__':main()
