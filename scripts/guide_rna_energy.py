"""Per-guide RNA energy tests with deterministic balanced controls.

The observed V-statistic uses pertpy 1.4.0 Edistance. The published unbiased
U-statistic is a positive affine transform for equal group sizes, so its
permutation ordering is identical. Batched arithmetic is checked against the
released implementation. Monte-Carlo p=(1+extreme)/(R+1), including ties.
No target-RNA decrease or response-selected input is used.
"""
from pathlib import Path
import argparse,csv,hashlib,json,os,socket,sys,time
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def atomic(path,value):
    path=Path(path);temp=path.with_suffix(path.suffix+'.receiving');temp.write_text(json.dumps(value,ensure_ascii=False,indent=2));temp.replace(path)


def energy_v(D,mask):
    mask=np.asarray(mask,dtype=bool);a=np.flatnonzero(mask);b=np.flatnonzero(~mask)
    return float(2*D[np.ix_(a,b)].mean()-D[np.ix_(a,a)].mean()-D[np.ix_(b,b)].mean())


def balanced_values(D,masks):
    """Same V statistic for balanced binary group masks, via batched sums."""
    n=D.shape[0]//2
    if D.shape!=(2*n,2*n) or not np.all(masks.sum(axis=1)==n):raise ValueError('Groups are not balanced')
    product=masks@D
    inside=np.einsum('ij,ij->i',product,masks)
    selected_row_sums=masks@D.sum(axis=1)
    between=selected_row_sums-inside
    return (4*between-D.sum(dtype=np.float64))/n**2


def unbiased_from_balanced_v(values,D):
    n=D.shape[0]//2
    if n<2 or D.shape!=(2*n,2*n):raise ValueError('Need equal groups of at least two')
    return ((2*n-1)/(2*(n-1)))*values-D.sum(dtype=np.float64)/(2*n*n*(n-1))


def monte_carlo_p(observed,null):
    # Treat numerical near-ties conservatively; exact identical data give p=1.
    tol=1e-10*max(1.,abs(float(observed)))
    extreme=int(np.count_nonzero(null>=observed-tol))
    return (extreme+1)/(len(null)+1),extreme


def test_pair(X,Y,seed,permutations,official=None,batch_size=256):
    from sklearn.metrics import pairwise_distances
    X=np.asarray(X,dtype=np.float64);Y=np.asarray(Y,dtype=np.float64)
    if X.shape!=Y.shape or len(X)<2:raise ValueError('Matched equal-N RNA groups required')
    if not np.isfinite(X).all() or not np.isfinite(Y).all():raise ValueError('Nonfinite RNA coordinates')
    n=len(X);D=pairwise_distances(np.vstack([X,Y]),metric='euclidean')
    if not np.allclose(D,D.T,atol=1e-12):raise ValueError('Pairwise matrix is not symmetric')
    mask=np.arange(2*n)<n
    observed=energy_v(D,mask) if official is None else float(official.from_precomputed(D,mask,total=D.sum(dtype=np.float64)))
    fast=float(balanced_values(D,mask[None,:].astype(float))[0])
    if not np.isclose(observed,fast,rtol=1e-9,atol=1e-9):raise ValueError('Observed statistic differs from released implementation')
    rng=np.random.default_rng(seed);template=mask.astype(np.float64);null=np.empty(permutations,dtype=np.float64)
    for start in range(0,permutations,batch_size):
        end=min(permutations,start+batch_size)
        masks=rng.permuted(np.broadcast_to(template,(end-start,2*n)),axis=1)
        null[start:end]=balanced_values(D,masks)
    p,extreme=monte_carlo_p(observed,null)
    null_u=unbiased_from_balanced_v(null,D)
    return dict(energy_biased=observed,energy_unbiased=float(unbiased_from_balanced_v(observed,D)),
                pvalue=p,extreme=extreme,permutations=permutations,n_each=n,
                null_u_mean=float(null_u.mean()),null_u_sd=float(null_u.std(ddof=1)),
                null_u_quantiles=np.quantile(null_u,[.025,.5,.975]).tolist()),null_u


def reference_rows(guide,guide_labels,roles,barcodes,n,seed):
    available=np.flatnonzero((roles=='non_targeting_control')&(guide_labels!=guide))
    if len(available)<n:raise ValueError('Insufficient independent NTC records')
    ranked=sorted(available,key=lambda i:hashlib.sha256((str(seed)+'|NTC_reference|'+str(barcodes[i])).encode()).digest())
    return np.asarray(ranked[:n],dtype=np.int64)


def calculate(cfg):
    import pertpy as pt
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    from importlib.metadata import version
    task=json.loads(Path(cfg['task_manifest']).read_text());op=Path(task['native_directory']);root=Path(task['canonical_output_root'])
    if os.getuid()!=2601 or socket.gethostname()!='rhesusbase':raise ValueError('Wrong execution owner/host')
    if sha(__file__)!=cfg['script_sha256']:raise ValueError('Fixed execution source differs')
    if version('pertpy')!='1.4.0':raise ValueError('Wrong released energy implementation')
    for p,h in cfg['module_SHA256'].items():
        if sha(p)!=h:raise ValueError('Library source differs')
    for p,h in cfg['input_SHA256'].items():
        if sha(p)!=h:raise ValueError('Frozen scientific input differs')
    root.mkdir(parents=True,exist_ok=True)
    lock=op/'active.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    st=Path('/proc/self/stat').read_text();ticks=int(st[st.rfind(')')+2:].split()[19]);os.write(fd,json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=ticks,argv=sys.argv,commit=cfg['commit'])).encode());os.close(fd)
    contract=digest({k:v for k,v in cfg.items() if k!='script_sha256'})
    def progress(stage,**kw):atomic(op/'status.json',dict(stage=stage,at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kw))
    try:
        progress('load_accepted_RNA_PC_checkpoint')
        with np.load(cfg['input_path'],allow_pickle=False) as z:
            pc=z['PCA'];barcodes=z['cell_barcode'];labels=z['guide'];roles=z['role'];RNA_ids=z['RNA_ids'];totals=z['total_RNA']
        if pc.shape!=(71664,40) or len(set(barcodes))!=71664:raise ValueError('RNA cohort identity differs')
        if hashlib.sha256(json.dumps(RNA_ids.tolist()).encode()).hexdigest()!=task['intended_RNA_inventory_SHA256']:raise ValueError('Intended ordered RNA inventory differs')
        cells=pd.read_csv(cfg['all_cells_path'],dtype={'cell_barcode':str,'guide_id':str});counts=cells.groupby('guide_id').size()
        expected=cells[cells.guide_id.map(counts).ge(10)]
        if expected.cell_barcode.tolist()!=barcodes.tolist() or expected.guide_id.tolist()!=labels.tolist():raise ValueError('Full eligible barcode/guide axes differ')
        catalogue=pd.read_csv(cfg['catalogue_path'],keep_default_na=False);catalogue=catalogue.set_index('guide',verify_integrity=True)
        groups=[]
        for g in np.unique(labels):
            rows=np.flatnonzero(labels==g);rr=np.unique(roles[rows])
            if len(rr)!=1 or str(rr[0])!=str(catalogue.loc[g,'source_control_type']):raise ValueError('Guide/design role differs')
            groups.append(dict(guide=str(g),role=str(rr[0]),rows=rows,n=len(rows)))
        groups.sort(key=lambda g:(g['role']!='non_targeting_control',hashlib.sha256(g['guide'].encode()).hexdigest()))
        if len(groups)!=785:raise ValueError('Incomplete guide inventory')
        official=pt.tl.Distance('edistance',obsm_key='X_pca',cell_wise_metric='euclidean').metric_fct
        (root/'guides').mkdir(exist_ok=True);results=[]
        progress('NTC_diagnostics_then_targeting',completed=0,total=785,permutations=cfg['permutations'])
        for i,g in enumerate(groups):
            file=root/'guides'/f'{i:04d}.json';null_path=file.with_suffix('.npz')
            ref=reference_rows(g['guide'],labels,roles,barcodes,g['n'],cfg['reference_seed'])
            sample_sha=digest(dict(query=barcodes[g['rows']].tolist(),reference=barcodes[ref].tolist()))
            if file.exists():
                row=json.loads(file.read_text())
                if row['contract_sha256']!=contract or row['sample_sha256']!=sample_sha or row['guide']!=g['guide']:raise ValueError('Checkpoint source/method differs')
                if not null_path.exists() or sha(null_path)!=row['null_SHA256']:raise ValueError('Completed null checkpoint differs')
            else:
                if null_path.exists():raise ValueError('Unaccepted partial guide artifact; preserve and recover explicitly')
                started=time.monotonic();seed=int.from_bytes(hashlib.sha256((str(cfg['permutation_seed'])+'|'+g['guide']).encode()).digest()[:8],'little')
                stat,null=test_pair(pc[g['rows']],pc[ref],seed,cfg['permutations'],official)
                temp=null_path.with_name(null_path.name+'.receiving')
                with temp.open('wb') as f:np.savez_compressed(f,null_energy_unbiased=null)
                temp.replace(null_path)
                row=dict(guide=g['guide'],role=g['role'],target_gene_name=str(catalogue.loc[g['guide'],'target_gene_name']),target_gene_id=str(catalogue.loc[g['guide'],'target_gene_id']),
                    **stat,seconds=time.monotonic()-started,contract_sha256=contract,sample_sha256=sample_sha,
                    query_barcodes=barcodes[g['rows']].tolist(),reference_barcodes=barcodes[ref].tolist(),
                    null_SHA256=sha(null_path),permutation_seed=seed,reference_seed=cfg['reference_seed'],median_RNA_UMI=float(np.median(totals[g['rows']])))
                atomic(file,row)
            results.append(row)
            ntc=[r for r in results if r['role']=='non_targeting_control']
            if len(ntc)==106 and (i==105):
                frac=sum(r['pvalue']<.05 for r in ntc)/106
                diagnostic=dict(NTC_guides=106,raw_p_lt05_count=sum(r['pvalue']<.05 for r in ntc),raw_p_lt05_fraction=frac,
                    review_required=frac>cfg['NTC_diagnostic_flag_fraction'],flag_is_diagnostic_not_formal_FDR_guarantee=True)
                atomic(root/'NTC_diagnostic.json',diagnostic)
                if diagnostic['review_required']:
                    progress('NTC_review_required_before_target_tests',completed=106,total=785,diagnostic=diagnostic)
                    return
            progress('guide_tests_running',completed=len(results),total=785,last_guide=g['guide'],last_seconds=row['seconds'])
        p=np.array([r['pvalue'] for r in results]);reject,adj,_,_=multipletests(p,alpha=.05,method='fdr_bh')
        small=[]
        for r,a,flag in zip(results,adj,reject):
            small.append({k:v for k,v in r.items() if k not in ['query_barcodes','reference_barcodes']})
            small[-1].update(pvalue_adj=float(a),significant_adj=bool(flag))
        frame=pd.DataFrame(small);frame.to_csv(root/'tested_guides.csv',index=False)
        all_table=catalogue.reset_index().merge(frame,on=['guide'],how='left',suffixes=('_design',''))
        all_table['n_QC_records']=all_table.guide.map(counts).fillna(0).astype(int)
        all_table['assessment_status']=np.where(all_table['n_QC_records'].eq(0),'no_QC_records',np.where(all_table['n_QC_records'].lt(10),'below10_not_tested','tested_RNA_state_not_functional_KO'))
        all_table.to_csv(root/'all_design_guides.csv',index=False)
        diagnostic=json.loads((root/'NTC_diagnostic.json').read_text())
        summary=dict(task_id=task['task_id'],commit=cfg['commit'],contract_sha256=contract,records=71664,tested_guides=len(frame),all_design_guides=len(catalogue),
            calibration=diagnostic,multiple_testing='BH across all785 tested guide groups including NTC/intergenic diagnostics',alpha=.05,
            targeting_significant=int(((frame.role.isin(['targeting','positive_control']))&frame.significant_adj).sum()),NTC_significant=int(((frame.role=='non_targeting_control')&frame.significant_adj).sum()),
            estimator='published equal-N unbiased energy U-statistic; affine-equivalent permutation ranking to released pertpy V statistic',
            cell_distance='Euclidean, frozen pre-response RNA PCA40',balanced_reference='all query records, same-number NTC reference, entire own NTC guide excluded',
            pvalue='(1+number of null statistics >= observed)/(R+1), numerical ties conservative',permutations=cfg['permutations'],
            no_target_RNA_filter=True,no_BAM_or_original_counts_read=True,inference_scope=task['input_scope'],functional_KO_efficiency_assessed=False)
        atomic(root/'summary.json',summary)
        atomic(root/'result_manifest.json',{n:dict(bytes=(root/n).stat().st_size,sha256=sha(root/n)) for n in ['tested_guides.csv','all_design_guides.csv','summary.json','NTC_diagnostic.json']})
        progress('guide_RNA_tests_complete_pending_acceptance',completed=785,total=785)
    except BaseException as e:progress('failed_pending_recovery',error=repr(e));raise
    finally:lock.unlink()


def benchmark(output):
    import pertpy as pt
    from sklearn.metrics import pairwise_distances
    official=pt.tl.Distance('edistance',obsm_key='X_pca').metric_fct;rng=np.random.default_rng(421)
    checks=[]
    for n in [10,40,100]:
        data=rng.normal(size=(2*n,40));D=pairwise_distances(data)
        masks=rng.permuted(np.broadcast_to((np.arange(2*n)<n).astype(float),(32,2*n)),axis=1)
        fast=balanced_values(D,masks)
        slow=np.array([official.from_precomputed(D,m.astype(bool),total=D.sum(dtype=np.float64)) for m in masks])
        error=float(np.max(np.abs(fast-slow)))
        if error>1e-9:raise ValueError('Batched released-statistic equivalence failed')
        checks.append(dict(n=n,max_abs_error=error))
    timing=[]
    for n in [30,100,250,500]:
        X=rng.normal(size=(n,40));Y=rng.normal(size=(n,40));start=time.monotonic();result,_=test_pair(X,Y,482,1000,official)
        timing.append(dict(n_each=n,permutations=1000,seconds=time.monotonic()-start))
    proof=dict(released_metric_equivalence=checks,timing=timing,scientific_data_used=False)
    atomic(output,proof);print(json.dumps(proof))


if __name__=='__main__':
    parser=argparse.ArgumentParser();s=parser.add_subparsers(dest='command',required=True)
    a=s.add_parser('calculate');a.add_argument('--config',type=Path,required=True)
    b=s.add_parser('benchmark');b.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='calculate':calculate(json.loads(args.config.read_text()))
    else:benchmark(args.output)
