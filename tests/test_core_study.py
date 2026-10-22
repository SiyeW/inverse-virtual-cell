"""Source-protection, exact receipt coverage and non-executing study setup."""
import copy
import hashlib
import json
from pathlib import Path
import socket
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import manage_core_study as study


def put(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) if not isinstance(value, str) else value, encoding='utf-8')
    return dict(path=str(path), sha256=study.digest(path))


def fixture(tmp_path):
    delivery = tmp_path/'protected'
    files, verified, libraries = {}, {}, []
    for sample, role in [('RNA','Gene Expression'), ('GUIDE','CRISPR Guide Capture')]:
        keys = []
        for read in ['R1', 'R2']:
            key = f'fastqs/{sample}_S1_L001_{read}_001.fastq.gz'
            path = delivery/key
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b'fixture fastq bytes')
            s = path.stat()
            identity = dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)
            proof = dict(sha256=study.digest(path), md5=hashlib.md5(path.read_bytes()).hexdigest(), identity=identity)
            files[key] = identity
            verified[key] = dict(source=proof, target=copy.deepcopy(proof), vendor_md5=proof['md5'])
            keys.append(key)
        libraries.append(dict(fastq_id=sample, feature_type=role, fastq_dir=str(delivery/'fastqs'), receipt_keys=keys))
    manifest = put(delivery/'manifest.json', dict(target=str(delivery), files=files,
        vendor_md5={k:v['vendor_md5'] for k,v in verified.items()}))
    ledger = put(delivery/'ledger.json', dict(manifest_sha256=manifest['sha256'], verified=verified))
    feature = put(delivery/'features.csv', 'id,name\na,a\n')
    config = dict(schema_version=1, study_id='study_fixture', execution_host=socket.gethostname(),
        code_commit='1'*40, study_root=str(tmp_path/'new_study'), protected_roots=[str(delivery)],
        source_manifest=manifest, source_ledger=ledger, feature_reference=feature,
        create_bam=False, capture_units_confirmed=False,
        chips=[dict(chip_id='chip07', expected_fastq_files=4, expected_fastq_bytes=sum(x['bytes'] for x in files.values()), libraries=libraries)],
        joint=dict(chip_ids=['chip07'], retain_chip_identity=True, experiment_relation='not_confirmed'),
        reuse={}, reference=dict(expected_name='custom_reference', path=None, content_receipt=None),
        runtime=dict(cellranger=None, version='10.1.0', work_root=None, reserved_workspace_bytes=None, filesystem_probe_accepted=False))
    return config


def save_config(tmp_path, config):
    p = tmp_path/'config.json'
    put(p, config)
    return p


def clean_revision(monkeypatch):
    monkeypatch.setattr(study.subprocess, 'check_output', lambda args, **kwargs: '1'*40 if 'rev-parse' in args else '')


def ready(tmp_path, monkeypatch, config):
    ref = tmp_path/'reference'
    files = {}
    for name in ['reference.json','fasta/genome.fa','genes/genes.gtf','star/Genome']:
        p = ref/name
        put(p, '{}')
        files[name] = dict(bytes=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns, sha256=study.digest(p))
    receipt = put(tmp_path/'reference_receipt.json', dict(complete=True, hash_algorithm='sha256', reference_root=str(ref), files=files))
    work = tmp_path/'workspace'; work.mkdir()
    config['capture_units_confirmed'] = True
    config['reference'].update(path=str(ref), content_receipt=receipt)
    config['runtime'].update(cellranger=sys.executable, work_root=str(work), reserved_workspace_bytes=1, filesystem_probe_accepted=True)
    monkeypatch.setattr(study.subprocess, 'run', lambda *args, **kwargs: subprocess.CompletedProcess([], 0, stdout='cellranger 10.1.0', stderr=''))
    return config


def test_check_reports_real_blockers_and_never_hashes_fastq_bodies(tmp_path, monkeypatch):
    config = fixture(tmp_path)
    original = study.digest
    def guarded(path):
        assert not str(path).endswith('.fastq.gz')
        return original(path)
    monkeypatch.setattr(study, 'digest', guarded)
    report = study.check(config)
    assert not report['ready_for_launch_preparation']
    assert report['input_files'] == 4 and report['fastq_body_bytes_read'] == 0
    assert set(report['blockers']) == {'experimental_capture_units_not_confirmed','cellranger_not_available',
        'custom_transcriptome_reference_missing','working_filesystem_not_reserved'}
    assert not Path(config['study_root']).exists()


@pytest.mark.parametrize('problem', ['host','bam','duplicate_chip','unsafe_chip','relative_root','protected_child','protected_parent','discard_origin','missing_joint_chip','unrecorded_joint_chip'])
def test_invalid_scope_stops_before_creating_any_namespace(tmp_path, problem):
    c = fixture(tmp_path)
    if problem == 'host': c['execution_host'] = 'another-host'
    elif problem == 'bam': c['create_bam'] = True
    elif problem == 'duplicate_chip': c['chips'].append(copy.deepcopy(c['chips'][0]))
    elif problem == 'unsafe_chip': c['chips'][0]['chip_id'] = '../chip'
    elif problem == 'relative_root': c['study_root'] = 'relative'
    elif problem == 'protected_child': c['study_root'] = str(tmp_path/'protected'/'new')
    elif problem == 'protected_parent': c['study_root'] = str(tmp_path)
    elif problem == 'discard_origin': c['joint']['retain_chip_identity'] = False
    elif problem == 'missing_joint_chip': c['joint']['chip_ids'] = []
    else: c['joint']['chip_ids'].append('chip01')
    with pytest.raises(ValueError): study.validate(c)


@pytest.mark.parametrize('problem', ['missing_pair','missing_whole_library_files','duplicate_file','wrong_sample','traversal','wrong_type','size_changed','mtime_changed','vendor_changed','source_sha_changed','source_identity_changed','metadata_sha_changed'])
def test_missing_or_inconsistent_proofs_stop_without_recopy(tmp_path, problem):
    c = fixture(tmp_path)
    library = c['chips'][0]['libraries'][0]
    ledger = json.loads(Path(c['source_ledger']['path']).read_text())
    key = library['receipt_keys'][0]
    target = Path(c['protected_roots'][0])/key
    if problem == 'missing_pair': library['receipt_keys'].pop()
    elif problem == 'missing_whole_library_files': c['chips'][0]['libraries'].pop()
    elif problem == 'duplicate_file': library['receipt_keys'].append(key)
    elif problem == 'wrong_sample': library['fastq_id'] = 'other'
    elif problem == 'traversal': library['receipt_keys'][0] = '../escape.fastq.gz'
    elif problem == 'wrong_type': library['feature_type'] = 'unknown'
    elif problem == 'size_changed': target.write_bytes(b'changed')
    elif problem == 'mtime_changed':
        import os
        os.utime(target, ns=(target.stat().st_atime_ns, target.stat().st_mtime_ns + 1_000_000_000))
    elif problem == 'vendor_changed': ledger['verified'][key]['vendor_md5'] = '0'*32
    elif problem == 'source_sha_changed': ledger['verified'][key]['source']['sha256'] = '0'*64
    elif problem == 'source_identity_changed': ledger['verified'][key]['source']['identity']['bytes'] += 1
    else: c['source_manifest']['sha256'] = '0'*64
    if problem in ['vendor_changed','source_sha_changed','source_identity_changed']:
        c['source_ledger'] = put(Path(c['source_ledger']['path']), ledger)
    with pytest.raises((ValueError, KeyError)): study.check(c)
    assert not Path(c['study_root']).exists()


def test_only_inode_change_preserves_content_receipt_only_on_explicit_actual_cifs(tmp_path, monkeypatch):
    c = fixture(tmp_path)
    ledger = json.loads(Path(c['source_ledger']['path']).read_text())
    key = c['chips'][0]['libraries'][0]['receipt_keys'][0]
    ledger['verified'][key]['target']['identity']['inode'] += 100
    c['source_ledger'] = put(Path(c['source_ledger']['path']), ledger)
    with pytest.raises(ValueError,match='actual-CIFS'): study.check(c)
    c['source_identity_policy'] = 'cifs-inode-audit-v1'
    monkeypatch.setattr(study,'filesystem_kind',lambda path:'ext4')
    with pytest.raises(ValueError,match='actual-CIFS'): study.check(c)
    monkeypatch.setattr(study,'filesystem_kind',lambda path:'cifs')
    report = study.check(c)
    assert report['inode_only_changes'] == 1 and report['fresh_fastq_hashes'] == 0
    assert report['inventory'][0]['content_evidence'].startswith('prior_complete')


def test_init_creates_traceable_layout_no_raw_copies_and_protects_existing_attempt(tmp_path, monkeypatch):
    c = fixture(tmp_path); clean_revision(monkeypatch)
    before = {p:study.digest(p) for p in Path(c['protected_roots'][0]).rglob('*') if p.is_file()}
    report = study.initialize(save_config(tmp_path,c))
    out = Path(c['study_root'])
    assert not report['scientific_analysis_started'] and not (out/'owner.lock').exists()
    assert (out/'runs/chip07').is_dir() and (out/'joint').is_dir()
    assert (out/'records/input_inventory.json').is_file()
    assert not list(out.rglob('*.fastq.gz')) and not list(out.rglob('*.bam'))
    assert study.load_initialized(out) == c
    assert before == {p:study.digest(p) for p in before}
    with pytest.raises(FileExistsError): study.initialize(tmp_path/'config.json')
    with pytest.raises(ValueError,match='Cannot render'): study.render(c)
    assert not (out/'configs/cellranger').exists()
    with (out/'study.json').open('a') as stream: stream.write(' ')
    with pytest.raises(ValueError,match='setup changed'): study.load_initialized(out)


def test_dirty_or_wrong_revision_refuses_init(tmp_path, monkeypatch):
    c = fixture(tmp_path)
    monkeypatch.setattr(study.subprocess,'check_output',lambda *args,**kwargs:'different')
    with pytest.raises(ValueError,match='fixed clean'): study.initialize(save_config(tmp_path,c))
    assert not Path(c['study_root']).exists()


def test_render_is_no_bam_and_preserves_multiple_sequencing_batches(tmp_path, monkeypatch):
    c = ready(tmp_path,monkeypatch,fixture(tmp_path))
    manifest = json.loads(Path(c['source_manifest']['path']).read_text())
    ledger = json.loads(Path(c['source_ledger']['path']).read_text())
    keys = []
    for read in ['R1','R2']:
        key = f'resequenced/RNA_EXTRA_S2_L002_{read}_001.fastq.gz'
        path = Path(c['protected_roots'][0])/key
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(b'resequencing fixture')
        s = path.stat()
        identity = dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)
        proof = dict(sha256=study.digest(path),md5=hashlib.md5(path.read_bytes()).hexdigest(),identity=identity)
        manifest['files'][key] = identity
        manifest['vendor_md5'][key] = proof['md5']
        ledger['verified'][key] = dict(source=proof,target=copy.deepcopy(proof),vendor_md5=proof['md5'])
        keys.append(key)
    c['chips'][0]['libraries'].append(dict(fastq_id='RNA_EXTRA',feature_type='Gene Expression',fastq_dir=str(path.parent),receipt_keys=keys))
    c['chips'][0]['expected_fastq_files'] += 2
    c['chips'][0]['expected_fastq_bytes'] += 2*len(b'resequencing fixture')
    c['source_manifest'] = put(Path(c['source_manifest']['path']),manifest)
    ledger['manifest_sha256'] = c['source_manifest']['sha256']
    c['source_ledger'] = put(Path(c['source_ledger']['path']),ledger)
    out = Path(c['study_root']); (out/'configs').mkdir(parents=True)
    result = study.render(c)
    text = (out/'configs/cellranger/chip07.csv').read_text()
    assert 'create-bam,false' in text and 'Gene Expression' in text and 'CRISPR Guide Capture' in text
    assert 'RNA_EXTRA' in text and str(path.parent) in text
    assert result == dict(configs_written=1,workloads_started=0)
    with pytest.raises(FileExistsError): study.render(c)
    c['chips'][0]['libraries'].pop()
    with pytest.raises(ValueError,match='coverage/bytes'): study.check(c)


def test_reference_tamper_and_workspace_shortfall_block_render(tmp_path,monkeypatch):
    c = ready(tmp_path,monkeypatch,fixture(tmp_path))
    assert study.check(c)['ready_for_launch_preparation']
    Path(c['reference']['path'],'fasta/genome.fa').write_text('different')
    c['runtime']['reserved_workspace_bytes'] = 10**30
    c['runtime']['filesystem_probe_accepted'] = False
    report = study.check(c)
    assert not report['ready_for_launch_preparation']
    assert 'workspace_below_recorded_reservation' in report['blockers']
    assert any(x.startswith('reference_file_identity_differs') for x in report['blockers'])


def test_oversized_metadata_is_rejected_before_hashing(tmp_path,monkeypatch):
    c = fixture(tmp_path)
    monkeypatch.setattr(study,'MAX_METADATA',1)
    monkeypatch.setattr(study,'digest',lambda path:pytest.fail('Must check metadata budget before hashing'))
    with pytest.raises(ValueError,match='budget'): study.checked_metadata(c['source_manifest'])
