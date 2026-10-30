"""Fixed-size sensitivity of frozen guide-level RNA energy comparisons.

Membership is selected by coverage only, before evaluating outcomes. Original
P values are re-corrected on exactly the same cohort. This comparison measures
sample-size sensitivity, not successful control calibration or functional KO.
"""
from pathlib import Path
import argparse,hashlib,json,os,socket,sys,time
import numpy as np
import pandas as pd
from guide_rna_energy import atomic,digest,sha,reference_rows,test_pair


def fixed_rows(rows,barcodes,guide,n,seed):
    rows=np.asarray(rows,dtype=np.int64)
    if len(rows)<n or n<2:raise ValueError('Insufficient fixed-N coverage')
    if len(set(barcodes[rows]))!=len(rows):raise ValueError('Repeated query barcode')
    return np.asarray(sorted(rows,key=lambda i:hashlib.sha256(f'{seed}|query|{guide}|{barcodes[i]}'.encode()).digest())[:n],dtype=np.int64)


def same_cohort_baseline(frame,guides):
    from statsmodels.stats.multitest import multipletests
    if not frame.guide.is_unique or len(set(guides))!=len(guides):raise ValueError('Repeated guide')
    out=frame.set_index('guide').loc[guides].reset_index().copy()
    reject,q,_,_=multipletests(out.pvalue,alpha=.05,method='fdr_bh')
    out['pvalue_adj']=q;out['significant_adj']=reject
    return out


def run(cfg):
    import pertpy as pt
    from statsmodels.stats.multitest import multipletests
    assert os.getuid()==2601 and socket.gethostname()=='rhesusbase'
    op=Path(cfg['native']);root=Path(cfg['output']);root.mkdir(exist_ok=True)
    for p,h in cfg['input_SHA256'].items():assert sha(p)==h,p
    base=json.loads(Path(cfg['original_config']).read_text())
    assert sha(Path(__file__).with_name('guide_rna_energy.py'))==base['script_sha256']
    for p,h in base['module_SHA256'].items():assert sha(p)==h,p
    lock=op/'n50.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    stat=Path('/proc/self/stat').read_text();ticks=int(stat[stat.rfind(')')+2:].split()[19]);os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,commit=cfg['commit'])).encode());os.close(fd)
    def progress(stage,**kw):atomic(op/'n50.status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**kw))
    try:
        with np.load(base['input_path'],allow_pickle=False) as z:pc=z['PCA'];bc=z['cell_barcode'];labels=z['guide'];roles=z['role']
        assert pc.shape==(71664,40) and len(set(bc))==71664
        original=pd.read_csv(cfg['baseline_table']);eligible=original.loc[original.n_each.ge(cfg['n'])].copy()
        eligible=eligible.sort_values('guide').reset_index(drop=True);guides=eligible.guide.tolist()
        assert len(guides)==cfg['expected_groups']
        controls=roles=='non_targeting_control';assert int(controls.sum())==19342
        baseline=same_cohort_baseline(original,guides);baseline.to_csv(root/'same_cohort_all_cells.csv',index=False)
        official=pt.tl.Distance('edistance',obsm_key='X_pca',cell_wise_metric='euclidean').metric_fct
        contract=digest(cfg);chunkdir=root/'blocks';chunkdir.mkdir(exist_ok=True);results=[]
        for start in range(0,len(guides),10):
            chunk=guides[start:start+10];p=chunkdir/f'{start:04d}.json'
            if p.exists():
                saved=json.loads(p.read_text());assert saved['contract_sha256']==contract and [r['guide'] for r in saved['rows']]==chunk
                results.extend(saved['rows']);progress('fixed_N_tests',completed=len(results),total=len(guides));continue
            rows=[]
            for guide in chunk:
                all_rows=np.flatnonzero(labels==guide);prior=eligible.loc[eligible.guide.eq(guide)].iloc[0];assert len(all_rows)==int(prior.n_each)
                query=fixed_rows(all_rows,bc,guide,cfg['n'],cfg['query_seed'])
                ref=reference_rows(guide,labels,roles,bc,cfg['n'],base['reference_seed'])
                assert not set(bc[query]).intersection(bc[ref]) and np.all(controls[ref]) and np.all(labels[ref]!=guide)
                seed=int.from_bytes(hashlib.sha256((str(base['permutation_seed'])+'|'+guide).encode()).digest()[:8],'little')
                statistic,_=test_pair(pc[query],pc[ref],seed,base['permutations'],official)
                rows.append(dict(guide=guide,role=prior.role,target_gene_name=prior.target_gene_name,original_n=int(prior.n_each),**statistic,
                    query_barcodes=bc[query].tolist(),reference_barcodes=bc[ref].tolist(),sample_SHA256=digest(dict(query=bc[query].tolist(),reference=bc[ref].tolist())),permutation_seed=seed))
            atomic(p,dict(contract_sha256=contract,rows=rows));results.extend(rows);progress('fixed_N_tests',completed=len(results),total=len(guides))
        full=pd.DataFrame([{k:v for k,v in r.items() if k not in ['query_barcodes','reference_barcodes']} for r in results]);assert full.guide.tolist()==guides
        reject,q,_,_=multipletests(full.pvalue,alpha=.05,method='fdr_bh');full['pvalue_adj']=q;full['significant_adj']=reject;full.to_csv(root/'fixed50.csv',index=False)
        summaries={}
        for name,frame in [('same_cohort_all_cells',baseline),('fixed50',full)]:
            values={}
            for role,mask in [('NTC',frame.role.eq('non_targeting_control')),('targeting',frame.role.isin(['targeting','positive_control'])),('intergenic',frame.role.eq('intergenic_control'))]:
                g=frame.loc[mask];values[role]=dict(guides=len(g),raw_P_lt05=int(g.pvalue.lt(.05).sum()),BH_P_lt05=int(g.pvalue_adj.lt(.05).sum()),energy_median=float(g.energy_unbiased.median()),energy_quantiles=g.energy_unbiased.quantile([.1,.25,.5,.75,.9]).tolist())
            summaries[name]=values
        summary=dict(commit=cfg['commit'],n=cfg['n'],query_seed=cfg['query_seed'],cohort_groups=len(full),same_guide_membership=True,all_cells_not_recomputed=True,BH_on_same_groups=True,
            comparisons=summaries,NTC_source_review_still_required=True,power_decreases_with_downsampling=True,lower_Pass_fraction_not_proof_problem_fixed=True,
            inferential_scope='sample-size sensitivity in merged records; not target-specific causality or functional KO',null_arrays_not_stored=True,contract_sha256=contract)
        atomic(root/'summary.json',summary);atomic(root/'manifest.json',{n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in ['same_cohort_all_cells.csv','fixed50.csv','summary.json']})
        progress('complete_pending_interpretation',completed=len(full),total=len(full))
    except BaseException as exc:progress('failed_pending_recovery',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);a=parser.parse_args();run(json.loads(a.config.read_text()))
