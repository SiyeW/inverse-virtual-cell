"""Finite NTC-only reference sensitivity; no targeting outcomes are examined."""
from pathlib import Path
import argparse,hashlib,json,os,time,sys
import numpy as np
from guide_rna_energy import test_pair,sha,atomic


def independent_reference(guide,labels,roles,barcodes,n,seed):
    rows=np.flatnonzero((roles=='non_targeting_control')&(labels!=guide))
    rows=sorted(rows,key=lambda i:hashlib.sha256((str(seed)+'|'+str(guide)+'|'+str(barcodes[i])).encode()).digest())
    if len(rows)<n:raise ValueError('Too few disjoint NTC records')
    return np.asarray(rows[:n],dtype=np.int64)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();cfg=json.loads(a.config.read_text())
    if os.getuid()!=2601:raise ValueError('Wrong owner')
    task=json.loads(Path(cfg['task_manifest']).read_text());op=Path(task['native_directory']);root=Path(task['canonical_output_root'])
    assert (op/'run.exit').read_text().strip()=='0' and not (op/'active.lock').exists()
    assert json.loads((op/'status.json').read_text())['stage']=='NTC_review_required_before_target_tests'
    if sha(cfg['input_path'])!=cfg['input_SHA256'][cfg['input_path']]:raise ValueError('RNA input differs')
    with np.load(cfg['input_path'],allow_pickle=False) as z:pc=z['PCA'];bc=z['cell_barcode'];labels=z['guide'];roles=z['role'];totalRNA=z['total_RNA']
    ntc_idx=np.flatnonzero(roles=='non_targeting_control');ntc=pc[ntc_idx].astype(float);mean=ntc.mean(axis=0);cov=np.cov(ntc,rowvar=False);inverse=np.linalg.inv(cov)
    groups=sorted(np.unique(labels[ntc_idx]),key=lambda g:hashlib.sha256(str(g).encode()).hexdigest());assert len(groups)==106
    import pertpy as pt
    official=pt.tl.Distance('edistance',obsm_key='X_pca').metric_fct
    result=[]
    cases=[('guide_specific_A',20261030),('guide_specific_B',20261101),('pooled_random_null',20261102)]
    destination=root/'reference_sensitivity';destination.mkdir(exist_ok=True)
    for case,seed in cases:
        file=destination/(case+'.json')
        if file.exists():result.append(json.loads(file.read_text()));continue
        rows=[]
        for i,g in enumerate(groups):
            query=np.flatnonzero(labels==g);n=len(query)
            if case=='pooled_random_null':
                rng=np.random.default_rng(int.from_bytes(hashlib.sha256((str(seed)+str(g)).encode()).digest()[:8],'little'))
                draw=rng.choice(ntc_idx,size=2*n,replace=False);query=draw[:n];reference=draw[n:]
            else:reference=independent_reference(g,labels,roles,bc,n,seed)
            perm_seed=int.from_bytes(hashlib.sha256((case+'|permutation|'+str(g)).encode()).digest()[:8],'little')
            stats,_=test_pair(pc[query],pc[reference],perm_seed,1999,official)
            delta=pc[reference].mean(axis=0)-mean
            rows.append(dict(guide=str(g),n=n,**{k:v for k,v in stats.items() if k in ['pvalue','energy_unbiased']},
                             query_median_RNA=float(np.median(totalRNA[query])),reference_median_RNA=float(np.median(totalRNA[reference])),reference_centroid_scaled_diagnostic=float(n*delta@inverse@delta)))
            atomic(op/'reference_diagnosis.status.json',dict(at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),case=case,completed=i+1,total=106))
        item=dict(case=case,seed=seed,permutations=1999,NTC_guides=106,raw_p_lt05_count=sum(r['pvalue']<.05 for r in rows),rows=rows)
        atomic(file,item);result.append(item)
    barcode_index={b:i for i,b in enumerate(bc)};old=[]
    for i in range(106):
        row=json.loads((root/'guides'/f'{i:04d}.json').read_text());ref=np.array([barcode_index[b] for b in row['reference_barcodes']]);delta=pc[ref].mean(axis=0)-mean
        old.append(float(len(ref)*delta@inverse@delta))
    summary=dict(cases=[dict(case=r['case'],raw_p_lt05_count=r['raw_p_lt05_count'],NTC_guides=106,pvalue_quantiles=np.quantile([s['pvalue'] for s in r['rows']],[0,.25,.5,.75,1]).tolist(),reference_centroid_scaled_diagnostic_quantiles=np.quantile([s['reference_centroid_scaled_diagnostic'] for s in r['rows']],[0,.25,.5,.75,1]).tolist()) for r in result],original_shared_reference_diagnostic_quantiles=np.quantile(old,[0,.25,.5,.75,1]).tolist(),diagnostic_not_new_response_filter=True,no_targeting_groups_examined=True)
    atomic(destination/'summary.json',summary);atomic(op/'reference_diagnosis.status.json',dict(at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),stage='complete',summary=summary));print(json.dumps(summary))
