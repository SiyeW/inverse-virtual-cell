from pathlib import Path
import gzip,json,sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from cache_chip01_htsjdk import verified_ready
from accept_imported_rna_cache import validate_manifest_inventory
from cache_guided_rna_alignments import sha_file,canonical_hash


def test_ready_cache_resume_rejects_payload_mutation(tmp_path):
    region=dict(region_id='R00000',chrom='chr1',start0=2,end0=30,guides=[])
    identity=dict(dataset='chip01',full_reference=[['chr1',100]])
    data=tmp_path/'R00000.jsonl.gz'
    with gzip.open(data,'wt') as f:f.write(json.dumps(dict(kind='RNA_SAM_cache_header',identity=identity,region=region))+'\n')
    ready=tmp_path/'R00000.ready.json'
    ready.write_text(json.dumps(dict(region=region,bytes=data.stat().st_size,sha256=sha_file(data))))
    assert verified_ready(ready,tmp_path,{'R00000':region},identity)['region']==region
    with data.open('ab') as f:f.write(b'changed')
    with pytest.raises(ValueError,match='content'):verified_ready(ready,tmp_path,{'R00000':region},identity)


def test_import_complete_inventory_does_not_accept_other_source_or_namespace():
    r=dict(region_id='R00000',chrom='chr1',start0=2,end0=30,guides=[])
    i=dict(dataset='chip01',reference=[['chr1',100]])
    receipt=dict(region=r,identity_sha256=canonical_hash(i),file='cache/R00000.jsonl.gz',bytes=100,sha256='a')
    m=dict(regions=[receipt],files={'cache/R00000.jsonl.gz':dict(bytes=100,sha256='a'),'source_header.sam':dict(bytes=10,sha256='b')})
    assert len(validate_manifest_inventory(m,{'R00000':r},i))==1
    with pytest.raises(ValueError,match='reference'):validate_manifest_inventory(m,{'R00000':r},dict(i,dataset='company'))
    m['files']['../source.bam']=dict(bytes=100,sha256='a')
    with pytest.raises(ValueError,match='roles'):validate_manifest_inventory(m,{'R00000':r},i)
