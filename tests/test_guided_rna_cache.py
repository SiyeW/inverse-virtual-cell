from pathlib import Path
from types import SimpleNamespace
import sys,json,gzip
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from cache_guided_rna_alignments import canonical_hash,sha_file,validate_region_receipt,samples_for_read


def test_full_reference_order_and_versions_are_material():
    a=dict(source='frozen',reference=[dict(name='chr1',length=100),dict(name='chr2',length=200)],code='v1')
    b=dict(a,reference=list(reversed(a['reference'])))
    assert canonical_hash(a)!=canonical_hash(b)
    assert canonical_hash(a)!=canonical_hash(dict(a,code='v2'))


def test_resume_requires_full_identity_not_matching_filename(tmp_path):
    identity=dict(task='t',reference=[['chr1',100]],assignment_sha='aa',code='v1')
    region=dict(region_id='R00000',chrom='chr1',start0=0,end0=80,guides=[])
    p=tmp_path/'R00000.jsonl.gz'
    with gzip.open(p,'wt') as f:f.write(json.dumps(dict(kind='RNA_SAM_cache_header',identity=identity,region=region))+'\n')
    item=dict(identity_sha256=canonical_hash(identity),region=region,bytes=p.stat().st_size,sha256=sha_file(p))
    validate_region_receipt(item,identity,region,p)
    with pytest.raises(ValueError,match='identity'):
        validate_region_receipt(item,dict(identity,assignment_sha='bb'),region,p)
    with p.open('ab') as f:f.write(b'broken')
    with pytest.raises(ValueError,match='content'):
        validate_region_receipt(item,identity,region,p)


def test_assignment_uses_exact_barcode_and_own_window():
    region=dict(guides=[dict(guide='A-1',sample_id='G0001',start0=200,end0=220)])
    def read(cb,start=195,end=225):
        return SimpleNamespace(has_tag=lambda k:k=='CB',get_tag=lambda k:cb,reference_start=start,reference_end=end)
    calls={'own':'A-1','other':'B-1','control':'NT-1'}
    assert samples_for_read(read('own'),region,calls,{'NT-1'})=='G0001'
    assert samples_for_read(read('control'),region,calls,{'NT-1'})=='NTC'
    assert samples_for_read(read('other'),region,calls,{'NT-1'}) is None
    assert samples_for_read(read('own',1000,1100),region,calls,{'NT-1'}) is None
    assert samples_for_read(read('own-1'),region,calls,{'NT-1'}) is None
