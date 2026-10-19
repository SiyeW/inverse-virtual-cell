import hashlib
import importlib.util
import json
import io
import threading
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location('delivery', Path(__file__).parents[1]/'scripts/copy_verified_delivery.py')
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


def setup_delivery(tmp_path):
    src,dst,control = (tmp_path/n for n in ('source','target','control'))
    src.mkdir(); control.mkdir()
    (src/'cellrangerresult').mkdir()
    (src/'cellrangerresult/a').write_bytes(b'abc'*100)
    m = dict(source=str(src),target=str(dst),control=str(control),cleanup_sources=False,
             files=d.tree_files(src),directories=['cellrangerresult'],vendor_md5={},
             total_files=1,total_bytes=300)
    p=control/'manifest.json'; p.write_text(json.dumps(m))
    return m,p,hashlib.sha256(p.read_bytes()).hexdigest()


def test_full_copy_and_resume_complete_keep_source(tmp_path):
    m,p,digest=setup_delivery(tmp_path)
    d.run(p,digest,m['control'])
    assert d.tree_files(m['source'])==m['files']
    assert (Path(m['target'])/'cellrangerresult/a').read_bytes()==b'abc'*100
    assert json.loads((Path(m['control'])/'final_acceptance.json').read_text())['complete']
    d.run(p,digest,m['control'])
    assert not (Path(m['control'])/'owner.lock').exists()


def test_partial_verified_before_resume(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest'
    src.write_bytes(b'abcdef')
    Path(str(dst)+'.receiving').write_bytes(b'abc')
    proof=d.copy_one(src,dst,d.identity(src),hashlib.md5(b'abcdef').hexdigest())
    assert dst.read_bytes()==b'abcdef'
    assert proof['source']['sha256']==proof['target']['sha256']


def test_identity_refreshes_cached_metadata(tmp_path,monkeypatch):
    p=tmp_path/'source';p.write_bytes(b'abc')
    real=Path.lstat
    class Cached:
        st_mode=real(p).st_mode
        st_size=3
        st_mtime_ns=1
        st_ino=0
        st_dev=0
    monkeypatch.setattr(Path,'lstat',lambda path:Cached() if path==p else real(path))
    refreshed=d.identity(p)
    assert refreshed['mtime_ns']==p.stat().st_mtime_ns
    assert refreshed['inode']==p.stat().st_ino


def test_atomic_json_recovers_from_temporary_permission_error(tmp_path,monkeypatch):
    p=tmp_path/'status.json';p.write_text('{"old": true}')
    original=d.os.replace
    calls=[]
    def transient(a,b):
        calls.append((a,b))
        if len(calls)<3:raise PermissionError(13,'temporary sharing denial')
        return original(a,b)
    monkeypatch.setattr(d.os,'replace',transient)
    monkeypatch.setattr(d.time,'sleep',lambda _:None)
    d.atomic_json(p,{'new':True})
    assert len(calls)==3 and json.loads(p.read_text())=={'new':True}


def test_atomic_json_persistent_permission_error_stops_and_preserves(tmp_path,monkeypatch):
    p=tmp_path/'status.json';p.write_text('{"old": true}')
    calls=[]
    def denied(*args):
        calls.append(args);raise PermissionError(13,'persistent denial')
    monkeypatch.setattr(d.os,'replace',denied)
    monkeypatch.setattr(d.time,'sleep',lambda _:None)
    with pytest.raises(PermissionError):d.atomic_json(p,{'new':True})
    assert len(calls)==6 and json.loads(p.read_text())=={'old':True}
    assert json.loads((tmp_path/'status.json.tmp').read_text())=={'new':True}


def test_destination_identity_remap_requires_fresh_hash(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc')
    frozen=d.identity(src);proof=d.copy_one(src,dst,frozen)
    proof['target']['identity']['inode']=-1
    calls=[]
    assert d.revalidate_published(src,dst,frozen,proof,lambda n,t:calls.append(n))
    assert calls[-1]==3 and proof['target']['identity']==d.identity(dst)
    assert len(proof['identity_revalidations'])==1


def test_revalidation_cannot_waive_same_size_content_change(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc')
    frozen=d.identity(src);proof=d.copy_one(src,dst,frozen)
    dst.write_bytes(b'xyz');proof['target']['identity']['inode']=-1
    with pytest.raises(RuntimeError,match='content changed'):
        d.revalidate_published(src,dst,frozen,proof)
    assert dst.read_bytes()==b'xyz' and src.read_bytes()==b'abc'


def test_bad_partial_never_appended(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest'
    src.write_bytes(b'abcdef')
    partial=Path(str(dst)+'.receiving');partial.write_bytes(b'bad')
    with pytest.raises(RuntimeError,match='prefix mismatch'):
        d.copy_one(src,dst,d.identity(src))
    assert partial.read_bytes()==b'bad' and not dst.exists()


def test_empty_partial_can_resume(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc')
    Path(str(dst)+'.receiving').write_bytes(b'')
    d.copy_one(src,dst,d.identity(src))
    assert dst.read_bytes()==b'abc'


def test_vendor_failure_retains_partial_not_final(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc')
    with pytest.raises(RuntimeError,match='Vendor MD5'):
        d.copy_one(src,dst,d.identity(src),'0'*32)
    assert src.read_bytes()==b'abc' and not dst.exists()
    assert Path(str(dst)+'.receiving').read_bytes()==b'abc'


def test_hash_failure_never_publishes(tmp_path,monkeypatch):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc')
    original=d.hashes
    def changed(*args,**kwargs):
        result=original(*args,**kwargs);result['sha256']='0'*64;return result
    monkeypatch.setattr(d,'hashes',changed)
    with pytest.raises(RuntimeError,match='SHA256 mismatch'):
        d.copy_one(src,dst,d.identity(src))
    assert not dst.exists()


def test_no_overwrite_and_source_change(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'dest';src.write_bytes(b'abc');dst.write_bytes(b'existing')
    with pytest.raises(FileExistsError): d.copy_one(src,dst,d.identity(src))
    assert dst.read_bytes()==b'existing'
    frozen=d.identity(src);src.write_bytes(b'other')
    with pytest.raises(RuntimeError,match='source changed'): d.copy_one(src,tmp_path/'other',frozen)


def test_second_owner_and_manifest_tamper(tmp_path):
    m,p,digest=setup_delivery(tmp_path)
    lock=Path(m['control'])/'owner.lock';lock.write_text('existing')
    with pytest.raises(FileExistsError): d.run(p,digest,m['control'])
    assert lock.read_text()=='existing'
    with pytest.raises(ValueError,match='SHA256'): d.run(p,'0'*64,m['control'])


def test_new_source_and_extra_target_rejected(tmp_path):
    m,p,digest=setup_delivery(tmp_path)
    Path(m['target']).mkdir();(Path(m['target'])/'unexpected').write_text('retain')
    with pytest.raises(RuntimeError,match='Unexpected destination'): d.run(p,digest,m['control'])
    assert (Path(m['target'])/'unexpected').read_text()=='retain'
    (Path(m['source'])/'new').write_text('x')
    with pytest.raises(RuntimeError,match='Source scope'): d.run(p,digest,m['control'])


def test_existing_destination_fresh_hash_not_size(tmp_path):
    m,p,digest=setup_delivery(tmp_path)
    dest=Path(m['target'])/'cellrangerresult/a';dest.parent.mkdir(parents=True)
    dest.write_bytes(b'xyz'*100)
    with pytest.raises(RuntimeError,match='mismatch'): d.run(p,digest,m['control'])
    assert dest.read_bytes()==b'xyz'*100


@pytest.mark.parametrize('rel', ['../escape','/absolute','x/../../escape',r'x\escape'])
def test_path_escape(tmp_path,rel):
    with pytest.raises(ValueError):d.safe_path(tmp_path,rel)


def test_symlink_rejected(tmp_path):
    if not hasattr(__import__('os'),'symlink'):pytest.skip('no symlink support')
    root=tmp_path/'root';root.mkdir();(root/'a').write_text('x')
    try:(root/'link').symlink_to(root/'a')
    except OSError:pytest.skip('symlink privilege unavailable')
    with pytest.raises(ValueError):d.safe_path(root,'link')


@pytest.mark.parametrize('parallel',[False,True])
def test_prefix_hash_order_and_unaligned_resume(tmp_path,monkeypatch,parallel):
    monkeypatch.setattr(d,'CHUNK',2)
    src,dst=tmp_path/'source',tmp_path/'target'
    src.write_bytes(b'abcdefghi');partial=Path(str(dst)+'.receiving');partial.write_bytes(b'abcde')
    events=[]
    proof=d.copy_one(src,dst,d.identity(src),hashlib.md5(b'abcdefghi').hexdigest(),
                     lambda *args:events.append(args),parallel_prefix=parallel)
    assert dst.read_bytes()==b'abcdefghi' and src.read_bytes()==b'abcdefghi'
    assert [(n,total) for stage,n,total in events if stage=='verify_resume_prefix']==[(2,5),(4,5),(5,5)]
    assert proof['source']['sha256']==proof['target']['sha256']==hashlib.sha256(b'abcdefghi').hexdigest()


def test_parallel_reads_really_overlap_and_hash_in_source_order(monkeypatch):
    monkeypatch.setattr(d,'CHUNK',2)
    rendezvous=threading.Barrier(2)
    class Paired(io.BytesIO):
        def read(self,n):
            assert n<=2
            rendezvous.wait(timeout=2)
            return super().read(n)
    sha,md5=hashlib.sha256(),hashlib.md5();events=[]
    d.verify_prefix(Paired(b'abcdef'),Paired(b'abcdef'),6,sha,md5,
                    lambda *a:events.append(a),True)
    assert sha.hexdigest()==hashlib.sha256(b'abcdef').hexdigest()
    assert md5.hexdigest()==hashlib.md5(b'abcdef').hexdigest()
    assert [a[1] for a in events]==[2,4,6]


@pytest.mark.parametrize('previous',[b'bad',b'ab'])
def test_parallel_prefix_mismatch_or_short_read_stops(previous):
    sha,md5=hashlib.sha256(),hashlib.md5()
    with pytest.raises(RuntimeError,match='prefix mismatch'):
        d.verify_prefix(io.BytesIO(b'abc'),io.BytesIO(previous),3,sha,md5,lambda *a:None,True)
    assert sha.hexdigest()==hashlib.sha256(b'').hexdigest()


@pytest.mark.parametrize('broken_source',[False,True])
def test_parallel_read_error_propagates(broken_source):
    class Broken(io.BytesIO):
        def read(self,n):raise OSError('simulated read failure')
    a,b=Broken(b'abc'),io.BytesIO(b'abc')
    if not broken_source:a,b=b,a
    with pytest.raises(OSError,match='read failure'):
        d.verify_prefix(a,b,3,hashlib.sha256(),hashlib.md5(),lambda *a:None,True)


def test_parallel_callback_failure_never_appends(tmp_path,monkeypatch):
    monkeypatch.setattr(d,'CHUNK',2)
    src,dst=tmp_path/'source',tmp_path/'target';src.write_bytes(b'abcdef')
    partial=Path(str(dst)+'.receiving');partial.write_bytes(b'abcd')
    def stop(*args):raise RuntimeError('stop before append')
    with pytest.raises(RuntimeError,match='stop before append'):
        d.copy_one(src,dst,d.identity(src),callback=stop,parallel_prefix=True)
    assert partial.read_bytes()==b'abcd' and not dst.exists()


def test_parallel_vendor_failure_retains_partial(tmp_path):
    src,dst=tmp_path/'source',tmp_path/'target';src.write_bytes(b'abcdef')
    partial=Path(str(dst)+'.receiving');partial.write_bytes(b'abc')
    with pytest.raises(RuntimeError,match='Vendor MD5'):
        d.copy_one(src,dst,d.identity(src),'0'*32,parallel_prefix=True)
    assert not dst.exists() and partial.read_bytes()==b'abcdef'


def test_parallel_run_keeps_lock_manifest_and_final_hash_gates(tmp_path):
    m,p,digest=setup_delivery(tmp_path)
    target=Path(m['target'])/'cellrangerresult/a';target.parent.mkdir(parents=True)
    Path(str(target)+'.receiving').write_bytes(b'abc'*30)
    d.run(p,digest,m['control'],parallel_prefix=True)
    final=json.loads((Path(m['control'])/'final_acceptance.json').read_text())
    ledger=json.loads((Path(m['control'])/'ledger.json').read_text())
    assert final['complete'] and final['parallel_prefix'] and not (Path(m['control'])/'owner.lock').exists()
    assert ledger['verified']['cellrangerresult/a']['source']['sha256']==ledger['verified']['cellrangerresult/a']['target']['sha256']
    assert d.tree_files(m['source'])==m['files']


def test_parallel_final_target_hash_failure_never_publishes(tmp_path,monkeypatch):
    src,dst=tmp_path/'source',tmp_path/'target';src.write_bytes(b'abcdef')
    partial=Path(str(dst)+'.receiving');partial.write_bytes(b'abc')
    original=d.hashes
    def mismatch(*args,**kwargs):
        proof=original(*args,**kwargs);proof['sha256']='0'*64;return proof
    monkeypatch.setattr(d,'hashes',mismatch)
    with pytest.raises(RuntimeError,match='SHA256 mismatch'):
        d.copy_one(src,dst,d.identity(src),parallel_prefix=True)
    assert not dst.exists() and partial.read_bytes()==b'abcdef'


def test_parallel_empty_prefix_has_no_reads_or_callbacks():
    class NoReads:
        def read(self,n):raise AssertionError('empty prefix must not read')
    d.verify_prefix(NoReads(),NoReads(),0,hashlib.sha256(),hashlib.md5(),
                    lambda *a:pytest.fail('unexpected callback'),True)


def test_parallel_prefix_never_reads_beyond_partial(monkeypatch):
    monkeypatch.setattr(d,'CHUNK',2)
    source,previous=io.BytesIO(b'abcdefghi'),io.BytesIO(b'abcde')
    d.verify_prefix(source,previous,5,hashlib.sha256(),hashlib.md5(),lambda *a:None,True)
    assert source.tell()==previous.tell()==5 and source.read()==b'fghi'
