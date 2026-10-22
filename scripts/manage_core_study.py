"""Create/check a source-preserving multi-chip study. Never starts a workload.

init --config private.json creates a new metadata/workflow namespace.
check --study PATH checks receipts/metadata and explicit execution prerequisites.
render-cellranger --study PATH writes configs only after check passes.
"""
import argparse
import csv
from datetime import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import socket
import subprocess


ROOT = Path(__file__).resolve().parents[1]
SHA = re.compile(r'[0-9a-f]{64}')
ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,95}')
READ = re.compile(r'_(S\d+)_(L\d{3})_([RI][12])_(\d{3})\.fastq\.gz$')
MAX_METADATA = 8 * 1024 * 1024


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_json(path):
    path = Path(path)
    if path.stat().st_size > MAX_METADATA:
        raise ValueError('Metadata file exceeds bounded inspection budget')
    return json.loads(path.read_text(encoding='utf-8'))


def write_new(path, value):
    """Exclusive creation: failures/evidence are never overwritten or resumed."""
    with Path(path).open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(value if isinstance(value, str) else json.dumps(value, ensure_ascii=False, indent=2) + '\n')
        stream.flush()
        os.fsync(stream.fileno())


def absolute(value):
    if not isinstance(value, str) or any(c in value for c in ['\x00', '\r', '\n']):
        raise ValueError('Invalid path')
    path = Path(value)
    if not path.is_absolute():
        raise ValueError('Absolute paths required on the execution host')
    return path.resolve()


def named(value):
    if not isinstance(value, str) or not ID.fullmatch(value):
        raise ValueError('Unsafe identifier')
    return value


def filesystem_kind(path):
    info = Path('/proc/self/mountinfo')
    if not info.is_file():
        return 'unknown'
    matches = []
    for line in info.read_text().splitlines():
        left, right = line.split(' - ', 1)
        mount = Path(left.split()[4].replace('\\040', ' '))
        if path == mount or path.is_relative_to(mount):
            matches.append((len(str(mount)), right.split()[0]))
    return max(matches)[1] if matches else 'unknown'


def file_spec(spec):
    if not isinstance(spec, dict) or not SHA.fullmatch(str(spec.get('sha256', ''))):
        raise ValueError('Frozen SHA256 file specification required')
    return absolute(spec['path'])


def checked_metadata(spec):
    path = file_spec(spec)
    if path.stat().st_size > MAX_METADATA:
        raise ValueError('Metadata file exceeds bounded inspection budget')
    if digest(path) != spec['sha256']:
        raise ValueError('Metadata SHA256 mismatch: ' + str(path))
    return read_json(path)


def validate(config):
    if config.get('schema_version') != 1:
        raise ValueError('Unsupported study schema')
    named(config['study_id'])
    if config.get('execution_host') != socket.gethostname():
        raise ValueError('Run on the recorded execution host, not a different workstation')
    if not re.fullmatch(r'[0-9a-f]{40}', config.get('code_commit', '')):
        raise ValueError('A fixed pushed Git commit is required')
    out = absolute(config['study_root'])
    protected = [absolute(p) for p in config['protected_roots']]
    if not protected or out == ROOT or out.is_relative_to(ROOT) or ROOT.is_relative_to(out):
        raise ValueError('A separate external study root is required')
    if any(out == p or out.is_relative_to(p) or p.is_relative_to(out) for p in protected):
        raise ValueError('Study root overlaps a protected source/result')
    specs = [config[k] for k in ['source_manifest', 'source_ledger', 'feature_reference']]
    for spec in specs:
        p = file_spec(spec)
        if not any(p.is_relative_to(root) for root in protected):
            raise ValueError('Input must belong to an explicit protected root')
    chips = config['chips']
    ids = [named(c['chip_id']) for c in chips]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError('Distinct chip identifiers required')
    if config.get('create_bam') is not False:
        raise ValueError('This study contract forbids new BAM output')
    if config.get('source_identity_policy', 'strict') not in ['strict', 'cifs-inode-audit-v1']:
        raise ValueError('Unknown source identity policy')
    seen = set()
    for chip in chips:
        if type(chip.get('expected_fastq_files')) is not int or chip['expected_fastq_files'] <= 0:
            raise ValueError('Expected per-chip FASTQ coverage required')
        if type(chip.get('expected_fastq_bytes')) is not int or chip['expected_fastq_bytes'] <= 0:
            raise ValueError('Expected per-chip FASTQ bytes required')
        roles = set()
        for library in chip['libraries']:
            named(library['fastq_id'])
            role = library['feature_type']
            if role not in ['Gene Expression', 'CRISPR Guide Capture']:
                raise ValueError('Unexpected library type')
            roles.add(role)
            absolute(library['fastq_dir'])
            if not library['receipt_keys']:
                raise ValueError('Empty library')
            for key in library['receipt_keys']:
                p = Path(key)
                if p.is_absolute() or '..' in p.parts or '\\' in key or key in seen:
                    raise ValueError('Unsafe or reused FASTQ receipt key')
                seen.add(key)
        if roles != {'Gene Expression', 'CRISPR Guide Capture'}:
            raise ValueError('RNA and CRISPR libraries must both be retained')
    joint = config['joint']
    joint_ids = joint['chip_ids']
    if len(joint_ids) != len(set(joint_ids)) or not set(ids).issubset(joint_ids):
        raise ValueError('Joint scope must retain every new chip exactly once')
    for chip_id in joint_ids:
        named(chip_id)
    reused = config.get('reuse', {})
    if set(joint_ids) - set(ids) != set(reused):
        raise ValueError('Every extra joint chip needs an explicit reused-input record')
    for spec in reused.values():
        p = file_spec(spec['source'])
        if not any(p.is_relative_to(root) for root in protected):
            raise ValueError('Reused input must be protected')
    if joint.get('retain_chip_identity') is not True:
        raise ValueError('Cell origin cannot be discarded')
    return config


def input_inventory(config):
    """Use complete accepted hashes; stat FASTQs, never reread their bodies."""
    manifest = checked_metadata(config['source_manifest'])
    ledger = checked_metadata(config['source_ledger'])
    if ledger['manifest_sha256'] != config['source_manifest']['sha256']:
        raise ValueError('Ledger does not bind the frozen source manifest')
    delivery = absolute(manifest['target'])
    if not any(delivery.is_relative_to(absolute(p)) for p in config['protected_roots']):
        raise ValueError('Delivery is not protected')
    rows = []
    for chip in config['chips']:
        for library in chip['libraries']:
            pairs = {}
            for key in library['receipt_keys']:
                frozen = manifest['files'][key]
                receipt = ledger['verified'][key]
                source, target = receipt['source'], receipt['target']
                if source['identity'] != frozen:
                    raise ValueError('Source receipt does not match frozen manifest: ' + key)
                sha = source['sha256']
                if not SHA.fullmatch(sha) or sha != target['sha256']:
                    raise ValueError('Missing/different source/target SHA256: ' + key)
                md5 = source['md5']
                if not re.fullmatch(r'[0-9a-f]{32}', md5) or md5 != target['md5']:
                    raise ValueError('Missing/different source/target MD5: ' + key)
                if receipt.get('vendor_md5') not in [None, md5]:
                    raise ValueError('Company MD5 differs: ' + key)
                vendor = manifest.get('vendor_md5', {}).get(key)
                if vendor is not None and (vendor != md5 or receipt.get('vendor_md5') != vendor):
                    raise ValueError('Frozen company MD5 proof missing/different: ' + key)
                if target['identity']['bytes'] != frozen['bytes']:
                    raise ValueError('Frozen size differs: ' + key)
                path = (delivery/key).resolve()
                if not path.is_relative_to(delivery):
                    raise ValueError('Resolved FASTQ escapes delivery')
                expected_parent = absolute(library['fastq_dir'])
                match = READ.search(path.name)
                if path.parent != expected_parent or not match or path.name[:match.start()] != library['fastq_id']:
                    raise ValueError('FASTQ library path/identity conflict: ' + key)
                sample, lane, read, chunk = match.groups()
                pairs.setdefault((sample, lane, chunk), set()).add(read)
                stat = path.stat()
                current = dict(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns, device=stat.st_dev, inode=stat.st_ino)
                old = target['identity']
                if any(current[k] != old[k] for k in ['bytes', 'mtime_ns', 'device']):
                    raise ValueError('Current source stable identity differs: ' + key)
                inode_changed = current['inode'] != old['inode']
                if inode_changed and (config.get('source_identity_policy') != 'cifs-inode-audit-v1' or filesystem_kind(path) != 'cifs'):
                    raise ValueError('Source inode differs outside explicit actual-CIFS audit policy: ' + key)
                rows.append(dict(chip_id=chip['chip_id'], feature_type=library['feature_type'], path=str(path),
                    receipt_key=key, sha256=sha, md5=md5, vendor_md5=receipt.get('vendor_md5'),
                    identity=current, inode_only_change=inode_changed,
                    content_evidence='prior_complete_source_target_hashes_plus_current_metadata', fresh_hash=False))
            if any(not {'R1', 'R2'}.issubset(reads) for reads in pairs.values()):
                raise ValueError('Missing paired RNA/guide read')
        selected = [r for r in rows if r['chip_id'] == chip['chip_id']]
        if len(selected) != chip['expected_fastq_files'] or sum(r['identity']['bytes'] for r in selected) != chip['expected_fastq_bytes']:
            raise ValueError('Per-chip expected FASTQ coverage/bytes differ')
    return rows


def check(config):
    validate(config)
    rows = input_inventory(config)
    feature = file_spec(config['feature_reference'])
    if feature.stat().st_size > MAX_METADATA:
        raise ValueError('Feature reference exceeds small-input budget')
    if digest(feature) != config['feature_reference']['sha256']:
        raise ValueError('Guide feature reference SHA256 mismatch')
    blockers = []
    runtime, reference = config['runtime'], config['reference']
    if config.get('capture_units_confirmed') is not True:
        blockers.append('experimental_capture_units_not_confirmed')
    executable = runtime.get('cellranger')
    if not executable or not absolute(executable).is_file():
        blockers.append('cellranger_not_available')
    else:
        result = subprocess.run([executable, '--version'], capture_output=True, text=True, timeout=10, check=True)
        if not re.fullmatch(r'\d+\.\d+\.\d+', runtime['version']) or not re.search(r'(?<![\d.])'+re.escape(runtime['version'])+r'(?![\d.])', result.stdout + result.stderr):
            blockers.append('cellranger_version_differs')
    reference_path = reference.get('path')
    if not reference_path or not absolute(reference_path).is_dir():
        blockers.append('custom_transcriptome_reference_missing')
    else:
        for name in ['reference.json', 'fasta/genome.fa', 'genes/genes.gtf', 'star/Genome']:
            if not (absolute(reference_path)/name).is_file():
                blockers.append('reference_component_missing:' + name)
        proof = reference.get('content_receipt')
        if not proof:
            blockers.append('reference_content_receipt_missing')
        else:
            record = checked_metadata(proof)
            if record.get('complete') is not True or record.get('reference_root') != str(absolute(reference_path)):
                blockers.append('reference_content_receipt_not_accepted')
            # No directory-wide rehash: a prior full receipt must include these exact files.
            for name in ['reference.json', 'fasta/genome.fa', 'genes/genes.gtf', 'star/Genome']:
                item = record.get('files', {}).get(name, {})
                p = absolute(reference_path)/name
                if not SHA.fullmatch(str(item.get('sha256', ''))) or not p.is_file():
                    blockers.append('reference_file_proof_missing:' + name)
                elif any(getattr(p.stat(), attr) != item.get(key) for attr, key in [('st_size','bytes'),('st_mtime_ns','mtime_ns')]):
                    blockers.append('reference_file_identity_differs:' + name)
            if record.get('hash_algorithm') != 'sha256' or not record.get('files'):
                blockers.append('reference_receipt_contract_incomplete')
            entries = list(record.get('files', {}).items())
            if len(entries) > 256:
                raise ValueError('Reference receipt exceeds bounded file budget')
            for name, item in entries:
                relative = Path(name)
                p = (absolute(reference_path)/relative).resolve()
                if relative.is_absolute() or '..' in relative.parts or not p.is_relative_to(absolute(reference_path)):
                    raise ValueError('Unsafe reference receipt path')
                if not SHA.fullmatch(str(item.get('sha256', ''))) or not p.is_file():
                    blockers.append('reference_file_proof_missing:' + name)
                elif p.stat().st_size != item.get('bytes') or p.stat().st_mtime_ns != item.get('mtime_ns'):
                    blockers.append('reference_file_identity_differs:' + name)
    work = runtime.get('work_root')
    reservation = runtime.get('reserved_workspace_bytes')
    workspace = None
    if not work or not absolute(work).is_dir():
        blockers.append('working_filesystem_not_reserved')
    else:
        work_path = absolute(work)
        if work_path == ROOT or work_path.is_relative_to(ROOT) or ROOT.is_relative_to(work_path):
            raise ValueError('Working directory overlaps code')
        if any(work_path == absolute(p) or work_path.is_relative_to(absolute(p)) or absolute(p).is_relative_to(work_path) for p in config['protected_roots']):
            raise ValueError('Working directory overlaps protected data')
        free = shutil.disk_usage(work_path).free
        workspace = dict(path=str(work_path), free_bytes=free, reserved_bytes=reservation)
        if type(reservation) is not int or reservation <= 0:
            blockers.append('explicit_workspace_budget_missing')
        elif free < reservation:
            blockers.append('workspace_below_recorded_reservation')
        if runtime.get('filesystem_probe_accepted') is not True:
            blockers.append('working_filesystem_probe_pending')
    return dict(schema_version=1, study_id=config['study_id'], checked_at=datetime.now().astimezone().isoformat(),
        ready_for_launch_preparation=not blockers, blockers=blockers, input_files=len(rows), input_bytes=sum(x['identity']['bytes'] for x in rows),
        fastq_body_bytes_read=0, fresh_fastq_hashes=0, metadata_receipt_policy=config.get('source_identity_policy', 'strict'),
        inode_only_changes=sum(x['inode_only_change'] for x in rows), workspace=workspace,
        joint_ready=False, joint_gates=['all_chip_core_results_accepted','shared_count_and_guide_contract','experimental_design_confirmed','chip_stratified_controls'],
        workloads_started=0, inventory=rows)


def initialize(config_path):
    config = validate(read_json(config_path))
    report = check(config)
    current = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
    if current != config['code_commit'] or dirty:
        raise ValueError('Initialize from the fixed clean code commit after tests/push')
    out = absolute(config['study_root'])
    if not out.parent.is_dir():
        raise ValueError('Study parent must already exist')
    out.mkdir()  # Refuse existing or failed namespaces; no automatic resume.
    write_new(out/'owner.lock', dict(pid=os.getpid(), hostname=socket.gethostname(), operation='metadata_initialization'))
    try:
        for name in ['records', 'configs', 'runs', 'joint', 'reports']:
            (out/name).mkdir()
        write_new(out/'study.json', config)
        write_new(out/'records/input_inventory.json', report.pop('inventory'))
        write_new(out/'records/preflight.initial.json', report)
        for chip in config['chips']:
            (out/'runs'/chip['chip_id']).mkdir()
            write_new(out/'configs'/(chip['chip_id']+'.cellranger.plan.json'),
                dict(chip_id=chip['chip_id'], libraries=chip['libraries'], create_bam=False,
                    expected_reference=reference_name(config), execution_ready=False))
        write_new(out/'records/workflow.json', dict(new_chips=[c['chip_id'] for c in config['chips']],
            per_chip=['cellranger_no_bam','accept_cellranger','prepare','core_smoke','core_formal','accept_report'],
            joint=config['joint'], status='planned_not_executed'))
        write_new(out/'README.md', '# '+config['study_id']+'\n\n'
            'Immutable deliveries stay outside this study. No raw data are copied here.\n'
            'study.json binds scope/code/inputs; records/ holds receipts and checks; configs/ holds parameters; '
            'runs/<chip>/ holds distinct attempts; joint/ holds origin-preserving combined analyses; '
            'reports/ is only for accepted reader deliverables.\n\n'
            'Initialization is not Cell Ranger, a scientific result, or final acceptance. '
            'Do not run drafts or overwrite a failed attempt. See records/preflight.initial.json.\n')
        files = {p.relative_to(out).as_posix():dict(bytes=p.stat().st_size, sha256=digest(p))
                 for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'owner.lock'}
        write_new(out/'records/setup_manifest.json', files)
    finally:
        (out/'owner.lock').unlink()
    return dict(study=str(out), scientific_analysis_started=False, setup_files=len(files), preflight=report)


def reference_name(config):
    return config['reference']['expected_name']


def load_initialized(study):
    """Check frozen setup bytes before reading/rerendering a study."""
    study = absolute(str(study))
    records = read_json(study/'records/setup_manifest.json')
    if 'study.json' not in records:
        raise ValueError('Study config is not bound by setup receipt')
    for name, proof in records.items():
        relative = Path(name)
        path = (study/relative).resolve()
        if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(study):
            raise ValueError('Unsafe setup receipt path')
        if not path.is_file() or path.stat().st_size != proof['bytes'] or digest(path) != proof['sha256']:
            raise ValueError('Frozen study setup changed: ' + name)
    return read_json(study/'study.json')


def render(config):
    report = check(config)
    if report['blockers']:
        raise ValueError('Cannot render executable configuration: '+', '.join(report['blockers']))
    out = absolute(config['study_root'])/'configs'/'cellranger'
    out.mkdir()
    for chip in config['chips']:
        stream = io.StringIO(newline='')
        writer = csv.writer(stream, lineterminator='\n')
        for row in [['[gene-expression]'],['reference',config['reference']['path']],['create-bam','false'],[],
                    ['[feature]'],['reference',config['feature_reference']['path']],[],
                    ['[libraries]'],['fastq_id','fastqs','feature_types']]:
            writer.writerow(row)
        for library in chip['libraries']:
            writer.writerow([library['fastq_id'],library['fastq_dir'],library['feature_type']])
        write_new(out/(chip['chip_id']+'.csv'), stream.getvalue())
    return dict(configs_written=len(config['chips']), workloads_started=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    init = commands.add_parser('init'); init.add_argument('--config', type=Path, required=True)
    for name in ['check', 'render-cellranger']:
        p = commands.add_parser(name); p.add_argument('--study', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'init':
        result = initialize(args.config)
    else:
        config = load_initialized(args.study)
        result = check(config) if args.command == 'check' else render(config)
        result.pop('inventory', None)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
