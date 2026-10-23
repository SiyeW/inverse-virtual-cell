import copy
from pathlib import Path
import sys
import json
from types import SimpleNamespace
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_cellranger_chip as runner
import manage_core_study as study

def config(tmp_path):
    work=tmp_path/'work';work.mkdir();(work/'outputs').mkdir();(work/'control').mkdir()
    root=tmp_path/'study';(root/'configs/cellranger').mkdir(parents=True)
    c=dict(study_root=str(root),reference=dict(path='/readonly/reference'),feature_reference=dict(path='/readonly/capture.csv'),
        runtime=dict(cellranger='/software/cellranger',work_root=str(work),localcores=8,localmem_gb=64,
                     reserved_workspace_bytes=140000000000,minimum_free_bytes=40000000000,rna_chemistry='SC5P-PE-v3',no_secondary=True),
        chips=[dict(chip_id='chip02',libraries=[dict(fastq_id='RNA',fastq_dir='/readonly/rna',feature_type='Gene Expression'),
                                             dict(fastq_id='CRISPR',fastq_dir='/readonly/crispr',feature_type='CRISPR Guide Capture')])])
    path=root/'configs/cellranger/chip02.csv';path.write_text(study.cellranger_csv(c,c['chips'][0]),newline='')
    return c,path

def test_explicit_new_run_preserves_inputs_no_bam_no_forced_cells(tmp_path):
    c,p=config(tmp_path);before=p.read_bytes();plan=runner.prepare_plan(c,'chip02','attempt1')
    assert plan['argv'][1]=='multi' and '--disable-ui' in plan['argv']
    assert '--localcores' in plan['argv'] and '--force-cells' not in plan['argv']
    assert b'create-bam,false\n' in before and b'chemistry,SC5P-PE-v3\n' in before
    assert b'no-secondary,true\n' in before
    assert p.read_bytes()==before and not plan['output'].exists() and not plan['control'].exists()
    plan['output'].mkdir()
    with pytest.raises(ValueError,match='unused'):runner.prepare_plan(c,'chip02','attempt1')

def test_config_tamper_and_unknown_chip_stop(tmp_path):
    c,p=config(tmp_path);p.write_text(p.read_text().replace('create-bam,false','create-bam,true'))
    with pytest.raises(ValueError,match='differs'):runner.prepare_plan(c,'chip02','attempt1')
    with pytest.raises(ValueError,match='recorded'):runner.prepare_plan(c,'chip99','attempt1')

@pytest.mark.parametrize('key,value',[('localcores',True),('localmem_gb',0),('reserved_workspace_bytes',-1),('minimum_free_bytes',None)])
def test_invalid_budgets_stop(tmp_path,key,value):
    c,_=config(tmp_path);c['runtime'][key]=value
    with pytest.raises(ValueError,match='budgets'):runner.limits(c['runtime'])

def test_guard_uses_free_floor_and_owned_bytes(tmp_path):
    c,_=config(tmp_path);p=runner.limits(c['runtime'])
    assert runner.budget_ok(50_000_000_000,100_000_000_000,p)
    assert not runner.budget_ok(30_000_000_000,100_000_000_000,p)
    assert not runner.budget_ok(50_000_000_000,150_000_000_000,p)
    d=tmp_path/'output';d.mkdir();(d/'a').write_bytes(b'abc')
    assert runner.owned_bytes(d)==3

def mock_ready(tmp_path, monkeypatch):
    c,p=config(tmp_path);c['code_commit']='a'*40
    monkeypatch.setattr(study,'load_initialized',lambda *_:c)
    monkeypatch.setattr(study,'check',lambda *_:dict(blockers=[],inventory=[]))
    monkeypatch.setattr(runner,'check_authority',lambda *_:{})
    monkeypatch.setattr(runner.subprocess,'check_output',lambda args,**kw:'a'*40 if 'rev-parse' in args else '')
    monkeypatch.setattr(runner,'available_memory',lambda:200*1024**3)
    monkeypatch.setattr(runner.os,'sched_getaffinity',lambda _:set(range(16)),raising=False)
    monkeypatch.setattr(runner,'process_identity',lambda:dict(pid=123,uid=456,startticks='789'))
    monkeypatch.setattr(runner.shutil,'disk_usage',lambda _:SimpleNamespace(free=500_000_000_000))
    return c,p

def test_one_explicit_run_stops_at_pending_acceptance(tmp_path,monkeypatch):
    c,p=mock_ready(tmp_path,monkeypatch);before=p.read_bytes();calls=[]
    class FakeProcess:
        pid=222;returncode=0
        def __init__(self,argv,**kw):
            calls.append((argv,kw));out=Path(argv[-1])/'outs';out.mkdir(parents=True);(out/'fixture.txt').write_text('fixture')
        def poll(self):return 0
    monkeypatch.setattr(runner.subprocess,'Popen',FakeProcess)
    result=runner.run(Path(c['study_root']),'chip02','attempt1')
    assert result['stage']=='cellranger_complete_pending_acceptance' and result['all_reports_complete'] is False
    assert len(calls)==1 and calls[0][1]['start_new_session'] is True
    assert not (Path(c['runtime']['work_root'])/'active_cellranger.lock').exists()
    assert p.read_bytes()==before
    with pytest.raises(ValueError,match='unused'):runner.run(Path(c['study_root']),'chip02','attempt1')

def test_existing_owner_is_not_replaced(tmp_path,monkeypatch):
    c,_=mock_ready(tmp_path,monkeypatch);lock=Path(c['runtime']['work_root'])/'active_cellranger.lock'
    lock.write_text('existing owner')
    monkeypatch.setattr(runner.subprocess,'Popen',lambda *a,**k:pytest.fail('Must not launch'))
    with pytest.raises(FileExistsError):runner.run(Path(c['study_root']),'chip02','attempt1')
    assert lock.read_text()=='existing owner'

def test_disk_guard_interrupts_owned_group_once_and_keeps_attempt(tmp_path,monkeypatch):
    c,_=mock_ready(tmp_path,monkeypatch);signals=[];launched=[]
    class FakeProcess:
        pid=222;returncode=None
        def __init__(self,*a,**kw):launched.append(self)
        def poll(self):return self.returncode
        def wait(self,timeout):self.returncode=-2;return -2
    def disk(_):return SimpleNamespace(free=500_000_000_000 if not launched else 30_000_000_000)
    monkeypatch.setattr(runner.shutil,'disk_usage',disk)
    monkeypatch.setattr(runner.subprocess,'Popen',FakeProcess)
    monkeypatch.setattr(runner.os,'killpg',lambda *args:signals.append(args),raising=False)
    with pytest.raises(RuntimeError,match='Disk guard'):runner.run(Path(c['study_root']),'chip02','attempt1')
    state=json.loads((Path(c['runtime']['work_root'])/'control/chip02_attempt1/status.json').read_text())
    assert state['stage']=='stopped_preserved' and len(launched)==1 and signals[0][0]==222
    assert not (Path(c['runtime']['work_root'])/'active_cellranger.lock').exists()
