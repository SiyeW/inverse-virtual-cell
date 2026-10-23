"""Compare complete official tiny fixture outputs, never production datasets."""
import argparse, hashlib, json
from pathlib import Path
import h5py
import numpy as np

def compare_matrix(a,b):
    checked=[]
    with h5py.File(a,'r') as left,h5py.File(b,'r') as right:
        names=[];other=[]
        left.visititems(lambda n,o:names.append(n) if isinstance(o,h5py.Dataset) else None)
        right.visititems(lambda n,o:other.append(n) if isinstance(o,h5py.Dataset) else None)
        if names!=other:raise ValueError('Matrix dataset coverage differs')
        for name in names:
            x,y=left[name][...],right[name][...]
            if x.shape!=y.shape or not np.array_equal(x,y):raise ValueError('Matrix values/order differ: '+name)
            checked.append(dict(dataset=name,shape=list(x.shape),equal=True))
    return checked

def inventory(root):
    files=[];total=0
    for p in sorted(root.rglob('*')):
        if not p.is_file():continue
        size=p.stat().st_size;total+=size
        if size>100_000_000 or total>200_000_000 or len(files)>100:
            raise ValueError('Not a bounded tiny fixture output')
        h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
        files.append(dict(path=str(p.relative_to(root)),bytes=size,sha256=h.hexdigest()))
    if any(x['path'].endswith('.bam') for x in files):raise ValueError('Unexpected BAM')
    return files

def verify(baseline,candidate):
    for root in [baseline,candidate]:
        status=json.loads((root/'status.json').read_text())
        if status.get('stage')!='engineering_fixture_passed' or status.get('scientific_result') is not False:
            raise ValueError('Only completed official engineering fixtures admitted')
        log=root/'tiny_no_bam/_log'
        if log.stat().st_size>2_000_000:raise ValueError('Unexpected fixture log size')
        if 'Could not move output files' in log.read_text():raise ValueError('Runtime reported output publication errors')
        final=root/'tiny_no_bam/SC_RNA_COUNTER_CS/fork0/_outs'
        if final.stat().st_size>100000:raise ValueError('Unexpected final output metadata size')
        json.loads(final.read_text())
    old=inventory(baseline/'tiny_no_bam/outs');new=inventory(candidate/'tiny_no_bam/outs')
    accepted=json.loads((baseline/'fixture_output_manifest.json').read_text())
    if old!=accepted:raise ValueError('Baseline differs from its accepted output hashes')
    comparisons={}
    for filename in ['raw_feature_bc_matrix.h5','filtered_feature_bc_matrix.h5']:
        comparisons[filename]=compare_matrix(baseline/'tiny_no_bam/outs'/filename,candidate/'tiny_no_bam/outs'/filename)
    if [x['path'] for x in old]!=[x['path'] for x in new]:raise ValueError('Output file coverage differs')
    for filename in ['metrics_summary.csv']:
        if (baseline/'tiny_no_bam/outs'/filename).read_bytes()!=(candidate/'tiny_no_bam/outs'/filename).read_bytes():
            raise ValueError('QC metrics differ: '+filename)
    return dict(complete=True,scientific_result=False,matrices=comparisons,metrics_equal=True,
                candidate_files=new,baseline_files=old,no_bam=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);a=p.parse_args()
    result=verify(a.baseline,a.candidate)
    with a.receipt.open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ['candidate_files','baseline_files','matrices']}))

if __name__=='__main__':main()
