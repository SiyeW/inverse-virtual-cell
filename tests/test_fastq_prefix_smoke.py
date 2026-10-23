import gzip
import hashlib
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import prepare_fastq_prefix_smoke as smoke

def fixture(tmp_path, mismatch=False, malformed=False):
    paths=[]
    for read in [1,2]:
        p=tmp_path/f'R{read}.gz'
        with gzip.open(p,'wb') as f:
            for i in range(5):
                name=f'@read{i + (1 if mismatch and read==2 else 0)}/{read}\n'
                f.write(name.encode()+b'ACGT\n+\n'+(b'II\n' if malformed else b'IIII\n'))
        paths.append(p)
    return paths

def test_paired_prefix_exact_and_source_unchanged(tmp_path):
    a,b=fixture(tmp_path);before=[p.read_bytes() for p in [a,b]]
    oa,ob=tmp_path/'out1.gz',tmp_path/'out2.gz'
    rows=smoke.copy_pair(a,b,oa,ob,3)
    assert [p.read_bytes() for p in [a,b]]==before
    for source,out,row in zip([a,b],[oa,ob],rows):
        with gzip.open(source,'rb') as f:expected=b''.join(f.readline() for _ in range(12))
        with gzip.open(out,'rb') as f:assert f.read()==expected
        assert hashlib.sha256(expected).hexdigest()==row['source_prefix_uncompressed_sha256']
    with pytest.raises(ValueError,match='prior output'):smoke.copy_pair(a,b,oa,ob,3)

@pytest.mark.parametrize('kwargs,pattern',[({'mismatch':True},'IDs differ'),({'malformed':True},'quality')])
def test_invalid_pairs_fail(tmp_path,kwargs,pattern):
    a,b=fixture(tmp_path,**kwargs)
    with pytest.raises(ValueError,match=pattern):smoke.copy_pair(a,b,tmp_path/'a.gz',tmp_path/'b.gz',3)

@pytest.mark.parametrize('count',[0,True,200001])
def test_unbounded_or_invalid_prefix_rejected(tmp_path,count):
    a,b=fixture(tmp_path)
    with pytest.raises(ValueError,match='1..200000'):smoke.copy_pair(a,b,tmp_path/'a.gz',tmp_path/'b.gz',count)

def test_short_source_and_overlap_rejected(tmp_path):
    a,b=fixture(tmp_path)
    with pytest.raises(ValueError,match='overlap'):smoke.copy_pair(a,b,a,tmp_path/'other',3)
    with pytest.raises(ValueError,match='ended'):smoke.copy_pair(a,b,tmp_path/'a.gz',tmp_path/'b.gz',6)
