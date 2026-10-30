"""Prespecified control experiment preserving shared-reference dependence.

Only NTC records are used. Real guide queries are frozen at the prior fixed50
draw. Random pseudoqueries preserve the reference policy, so failure caused
by a particular shared reference is not mislabeled guide-specific biology.
"""
from pathlib import Path
import argparse,json,hashlib,os,socket,sys,time
import numpy as np
import pandas as pd
from guide_rna_energy import atomic,digest,sha,reference_rows,test_pair
from compare_energy_fixed_n import fixed_rows
from diagnose_energy_reference import independent_reference


def random_query(guide,labels,roles,barcodes,reference,n,seed):
    available=np.flatnonzero((roles=='non_targeting_control')&(labels!=guide))
    available=available[~np.isin(available,reference)]
    if len(available)<n:raise ValueError('Too few disjoint random controls')
    return np.asarray(sorted(available,key=lambda i:hashlib.sha256(f'{seed}|query_null|{guide}|{barcodes[i]}'.encode()).digest())[:n],dtype=np.int64)


def run(cfg):
    import pertpy as pt
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    op=Path(cfg['native']);root=Path(cfg['output']);root.mkdir(exist_ok=True)
    old=json.loads(Path(cfg['original_config']).read_text());n50=json.loads(Path(cfg['fixed50_config']).read_text())
    for p,h in cfg['input_SHA256'].items():assert sha(p)==h,p
    assert sha(Path(__file__).with_name('guide_rna_energy.py'))==old['script_sha256']
    for p,h in old['module_SHA256'].items():assert sha(p)==h,p
    fd=os.open(op/'control.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY);s=Path('/proc/self/stat').read_text();ticks=int(s[s.rfind(')')+2:].split()[19]);os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,commit=cfg['commit'])).encode());os.close(fd)
    def progress(stage,**kw):atomic(op/'control.status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**kw))
    try:
        with np.load(old['input_path'],allow_pickle=False) as z:pc=z['PCA'];bc=z['cell_barcode'];labels=z['guide'];roles=z['role']
        baseline=pd.read_csv(cfg['fixed50_table']);nc=baseline.loc[baseline.role.eq('non_targeting_control')].sort_values('guide');assert len(nc)==99
        guide_list=nc.guide.tolist();contract=digest(cfg)
        official=pt.tl.Distance('edistance',obsm_key='X_pca',cell_wise_metric='euclidean').metric_fct
        allrows=[]
        for case in cfg['cases']:
            name=case['name'];dest=root/name;dest.mkdir(exist_ok=True);rows=[]
            for start in range(0,99,10):
                file=dest/f'{start:03d}.json';guides=guide_list[start:start+10]
                if file.exists():
                    saved=json.loads(file.read_text());assert saved['contract']==contract and [x['guide'] for x in saved['rows']]==guides;rows.extend(saved['rows']);continue
                block=[]
                for guide in guides:
                    q=fixed_rows(np.flatnonzero(labels==guide),bc,guide,50,n50['query_seed'])
                    if case['reference']=='shared':ref=reference_rows(guide,labels,roles,bc,50,old['reference_seed'])
                    else:ref=independent_reference(guide,labels,roles,bc,50,case['reference_seed'])
                    if case['query']=='pooled_random':q=random_query(guide,labels,roles,bc,ref,50,case['query_seed'])
                    assert len(q)==len(ref)==50 and not set(q).intersection(ref) and np.all(roles[q]=='non_targeting_control') and np.all(roles[ref]=='non_targeting_control')
                    assert np.all(labels[ref]!=guide)
                    seed=int.from_bytes(hashlib.sha256((str(old['permutation_seed'])+'|'+guide).encode()).digest()[:8],'little')
                    stats,_=test_pair(pc[q],pc[ref],seed,old['permutations'],official)
                    block.append(dict(case=name,guide=guide,**stats,sample_sha=digest(dict(query=bc[q].tolist(),reference=bc[ref].tolist()))))
                atomic(file,dict(contract=contract,rows=block));rows.extend(block);progress('NTC_control_design',case=name,completed=len(rows),total=99)
            allrows.extend(rows)
        frame=pd.DataFrame(allrows);frame.to_csv(root/'NTC_control_design.csv',index=False)
        summary={}
        for name,g in frame.groupby('case',sort=False):
            assert len(g)==99;summary[name]=dict(tests=99,raw_P_lt05=int(g.pvalue.lt(.05).sum()),energy_median=float(g.energy_unbiased.median()),pvalue_quantiles=g.pvalue.quantile([.1,.25,.5,.75,.9]).tolist())
        summary['original_shared_real_query']=dict(tests=99,raw_P_lt05=int(nc.pvalue.lt(.05).sum()),energy_median=float(nc.energy_unbiased.median()))
        proof=dict(commit=cfg['commit'],cases=summary,all_cases_fixed50=True,query_real_draw_unchanged=True,original_results_preserved=True,targeting_rows_analyzed=False,
            raw_P_counts_dependent_not_binomial_FDR=True,contract_sha256=contract)
        atomic(root/'summary.json',proof);atomic(root/'manifest.json',{n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in ['NTC_control_design.csv','summary.json']});progress('complete_pending_interpretation',cases=len(cfg['cases']),tests=len(frame))
    except BaseException as e:progress('failed_pending_recovery',error=repr(e));raise
    finally:(op/'control.lock').unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);a=p.parse_args();run(json.loads(a.config.read_text()))
