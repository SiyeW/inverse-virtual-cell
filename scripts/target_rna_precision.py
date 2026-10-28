"""Conditional mean uncertainty and equal-n NTC variability for target RNA.

Reuses the accepted target-count cache and freezes every original point.
Intervals describe cell sampling within this observed cohort, not independent
experiments, DNA editing, calibrated DE tests, or a perturbation success rate.
"""
from pathlib import Path
import argparse,hashlib,io,json,os,signal,time,warnings
import numpy as np
import pandas as pd
import scipy
from scipy import sparse,stats
from assess_target_rna import public_inputs
from compare_agreement_overlay import sha,atomic
from perturbseq_core_effects import stable_rng
from public_count_cache import check_arrays,axis_hash


def guide_interval(values,reference,seed,resamples=2000):
    values=np.asarray(values,float);n=len(values)
    if n<2 or reference<=0 or not np.isfinite(values).all() or (values<0).any():raise ValueError('Invalid guide interval input')
    point=float(values.mean()/reference);detected=int((values>0).sum())
    result=dict(RNA_remaining_ratio=point,guide_detected_cells=detected,guide_ci_low=np.nan,guide_ci_high=np.nan,
                guide_ci_method='SciPy BCa cell bootstrap; original NTC mean fixed',guide_ci_reference_fixed=True,
                zero_detection_fraction_upper95=1-.05**(1/n) if detected==0 else np.nan)
    if detected==0:result['guide_ci_status']='all_zero_no_RNA_mean_interval';return result
    if np.ptp(values)==0:result['guide_ci_status']='constant_observations_no_interval';return result
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit=stats.bootstrap((values,),np.mean,n_resamples=resamples,batch=128,vectorized=True,
                            confidence_level=.95,method='BCa',rng=np.random.default_rng(seed))
    low,high=fit.confidence_interval
    if not np.isfinite([low,high]).all() or any(issubclass(w.category,stats.DegenerateDataWarning) for w in caught):
        result['guide_ci_status']='bootstrap_degenerate_no_interval';return result
    result.update(guide_ci_low=float(low/reference),guide_ci_high=float(high/reference),guide_ci_status='estimated')
    return result


def matched_ntc_range(values,n,seed,resamples=2000):
    """No replacement within each NTC group; its reference excludes that group."""
    values=np.asarray(values,float);m=len(values)
    if n<1 or m-n<2 or not np.isfinite(values).all() or (values<0).any() or values.sum()<=0:
        raise ValueError('Insufficient NTC cells or missing RNA reference')
    total=float(values.sum());rng=np.random.default_rng(seed);ratios=np.empty(resamples)
    zero_reference=0
    for j in range(resamples):
        use=rng.choice(m,n,replace=False,shuffle=False);s=float(values[use].sum());other=(total-s)/(m-n)
        if other<=0:ratios[j]=np.inf;zero_reference+=1
        else:ratios[j]=(s/n)/other
    low,median,high=np.quantile(ratios,[.025,.5,.975],method='inverted_cdf')
    return dict(NTC_equal_n_low=float(low),NTC_equal_n_median=float(median),NTC_equal_n_high=float(high),
                NTC_null_draws=resamples,NTC_null_zero_reference_draws=zero_reference,NTC_null_sample_size=n,
                NTC_null_reference_excludes_selected_cells=True,NTC_constant_observations=bool(np.ptp(values)==0),
                NTC_mean_relative_SE=float(values.std(ddof=1)/np.sqrt(m)/values.mean()))


def saved_target_block(source,entry,columns,cells,remap,cache,contract):
    """Persist exact selected counts before intervals; a resume needs no old payload."""
    cache=Path(cache);receipt=cache.with_suffix('.json');columns=np.asarray(columns,int)
    if receipt.exists():
        saved=json.loads(receipt.read_text())
        if saved['contract']!=contract or saved['source_sha256']!=entry['sha256'] or sha(cache)!=saved['sha256']:
            raise ValueError('Selected target block changed')
        with np.load(cache,allow_pickle=False) as z:
            matrix=sparse.csc_matrix((z['data'],z['rows'],z['ptr']),shape=tuple(z['shape']))
            if not np.array_equal(z['gene_columns'],columns):raise ValueError('Selected RNA axis changed')
        if matrix.shape!=(cells,len(columns)):raise ValueError('Selected cell axis changed')
        return matrix,False
    source=Path(source);before=source.stat();payload=source.read_bytes()
    if len(payload)!=entry['bytes'] or hashlib.sha256(payload).hexdigest()!=entry['sha256']:
        raise ValueError('Accepted raw count block changed')
    with np.load(io.BytesIO(payload),allow_pickle=False) as z:data,rows,ptr=z['data'],z['rows'],z['ptr']
    check_arrays(data,rows,ptr,len(remap));del payload
    vals=[];idx=[];pointers=[0]
    for col in columns:
        j=col-entry['first'];a,b=int(ptr[j]),int(ptr[j+1]);rr=remap[rows[a:b]];keep=rr>=0
        vals.append(data[a:b][keep]);idx.append(rr[keep]);pointers.append(pointers[-1]+int(keep.sum()))
    matrix=sparse.csc_matrix((np.concatenate(vals),np.concatenate(idx),pointers),shape=(cells,len(columns)))
    if (source.stat().st_size,source.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):raise ValueError('Source changed during read')
    tmp=cache.with_suffix('.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,data=matrix.data,rows=matrix.indices,ptr=matrix.indptr,shape=np.array(matrix.shape),gene_columns=columns)
    tmp.replace(cache)
    atomic(receipt,dict(contract=contract,source_sha256=entry['sha256'],sha256=sha(cache),gene_columns=columns.tolist()))
    return matrix,True


def check_original_point(values,ref,observed):
    use=np.asarray(values,float)
    if len(use)!=int(observed.cells) or int((use>0).sum())!=int(observed.detected_cells):raise ValueError('Original guide cells/detection differ')
    if not np.isclose(use.mean(),float(observed.mean_per_10000),atol=1e-11,rtol=1e-11):raise ValueError('Original guide mean differs')
    if not np.isclose(ref,float(observed.NTC_mean_per_10000),atol=1e-11,rtol=1e-11):raise ValueError('Original NTC mean differs')


def contracts(d,config):
    source=hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest()
    method=dict(source=source,seed=config['seed'],resamples=config['resamples'],method='BCa_fixed_reference_plus_disjoint_NTC_equal_n',scipy=scipy.__version__,commit=config['expected_commit'])
    return source,hashlib.sha256(json.dumps(method,sort_keys=True).encode()).hexdigest()


def dataset(d,config,out,progress):
    out=Path(out);out.mkdir(exist_ok=True);(out/'parts').mkdir(exist_ok=True);(out/'selected_counts').mkdir(exist_ok=True)
    baseline=Path(d['target_result']);point_path=baseline/'guide_target_RNA.csv'
    if sha(point_path)!=d['point_table_sha256']:raise ValueError('Original point table changed')
    original=pd.read_csv(point_path)
    for filename,expected in d.get('source_metadata_sha256',{}).items():
        if sha(Path(filename))!=expected:raise ValueError('Frozen source metadata changed')
    source_contract,contract=contracts(d,config)
    proof={};consumed_bytes=0
    if d['kind']=='public_cache':
        cells,genes,_,manifest,all_obs=public_inputs(d['original_dataset'])
        columns=pd.Index(genes.gene_id.astype(str)).get_indexer(original.target_gene_id.unique())
        if (columns<0).any():raise ValueError('Original target RNA absent')
        columns=np.sort(columns)
        remap=np.full(len(all_obs),-1,np.int32);remap[cells.source_row.to_numpy(int)]=np.arange(len(cells))
        entries=manifest['blocks']
    else:
        if sha(baseline/'target_counts.json')!=d['target_counts_receipt_sha256']:raise ValueError('Lab target receipt changed')
        receipt=json.loads((baseline/'target_counts.json').read_text())
        if sha(baseline/'target_counts.npz')!=d['target_counts_sha256'] or receipt['sha256']!=d['target_counts_sha256']:
            raise ValueError('Saved lab target counts changed')
        for name,expected in d['cache_metadata_sha256'].items():
            if sha(baseline/name)!=expected:raise ValueError('Saved lab cell/RNA metadata changed')
        cells=pd.read_csv(baseline/'cells.csv.gz',keep_default_na=False);genes=pd.read_csv(baseline/'genes.csv.gz',keep_default_na=False)
        panel=sparse.load_npz(baseline/'target_counts.npz');columns=np.array(receipt['positions'],int)
        if panel.shape!=(len(cells),len(columns)):raise ValueError('Lab cached axes differ')
        entries=[dict(first=int(columns[i]),end=int(columns[min(i+32,len(columns))-1])+1,file='lab_%05d'%i) for i in range(0,len(columns),32)]
        proof['lab_target_cache_sha256']=d['target_counts_sha256'];proof['raw_count_matrix_reopened']=False
    if len(cells)!=d['expected_cells'] or int(cells.role.eq('non_targeting_control').sum())!=d['expected_NTC_cells']:
        raise ValueError('Frozen cohort differs')
    if not cells.cell_barcode.is_unique or (cells.total_counts<=0).any():raise ValueError('Invalid cell identity/totals')
    factors=10000./cells.total_counts.to_numpy(float);ntc=cells.role.eq('non_targeting_control').to_numpy()
    group_rows=cells.groupby('guide_id',sort=False).indices
    parts=[];selected_genes=0
    gene_to_column=dict(zip(genes.gene_id.astype(str),range(len(genes))))
    wanted=set(original.loc[original.cells.ge(10)&original.NTC_mean_per_10000.gt(0),'target_gene_id'])
    for entry in entries:
        relevant=np.array([j for j in columns if entry['first']<=j<entry['end'] and str(genes.iloc[j].gene_id) in wanted],int)
        if not len(relevant):continue
        part=out/'parts'/('%05d.csv'%entry['first']);pr=part.with_suffix('.json')
        if pr.exists():
            r=json.loads(pr.read_text())
            if r['contract']!=contract or sha(part)!=r['sha256']:raise ValueError('Completed uncertainty part changed')
            parts.append(pd.read_csv(part));selected_genes+=len(relevant);continue
        if d['kind']=='public_cache':
            source=Path(d['original_dataset']['baseline'])/'raw_count_cache'/entry['file']
            matrix,did_read=saved_target_block(source,entry,relevant,len(cells),remap,out/'selected_counts'/('%05d.npz'%entry['first']),source_contract)
            if did_read:consumed_bytes+=entry['bytes']
        else:matrix=panel[:,np.searchsorted(columns,relevant)]
        records=[]
        for k,j in enumerate(relevant):
            gene=str(genes.iloc[j].gene_id);values=matrix[:,k].toarray().ravel().astype(float)*factors
            nv=values[ntc];reference=float(nv.mean());null_cache={}
            targets=original.loc[original.target_gene_id.eq(gene)&original.cells.ge(10)&original.NTC_mean_per_10000.gt(0)]
            for row in targets.itertuples(index=False):
                v=values[group_rows[row.guide]];check_original_point(v,reference,row)
                key=int(row.cells)
                if key not in null_cache:
                    seed=int(stable_rng(config['seed'],d['id']+gene+str(key)+'NTC').integers(0,2**32-1))
                    null_cache[key]=matched_ntc_range(nv,key,seed,config['resamples'])
                seed=int(stable_rng(config['seed'],d['id']+row.guide+'CI').integers(0,2**32-1))
                ci=guide_interval(v,reference,seed,config['resamples'])
                if not np.isclose(ci['RNA_remaining_ratio'],row.RNA_remaining_ratio,atol=1e-11,rtol=1e-11):raise ValueError('Original ratio changed')
                records.append(dict(guide=row.guide,target_gene_id=gene,**ci,**null_cache[key]))
            selected_genes+=1
        frame=pd.DataFrame(records);frame.to_csv(part,index=False)
        atomic(pr,dict(contract=contract,sha256=sha(part),gene_columns=relevant.tolist(),source_entry=entry))
        parts.append(frame);progress('interval_blocks',completed_target_RNA=selected_genes,new_public_cache_bytes=consumed_bytes)
    precision=pd.concat(parts,ignore_index=True)
    if precision.guide.duplicated().any():raise ValueError('Duplicated interval guide')
    # Original point columns are retained exactly, not recalculated for plotting.
    precision=precision.drop(columns=['target_gene_id','RNA_remaining_ratio'])
    result=original.merge(precision,on='guide',how='left',validate='one_to_one')
    eligible=original.cells.ge(10)&original.NTC_mean_per_10000.gt(0)
    if result.loc[eligible,'guide_ci_status'].isna().any():raise ValueError('Missing uncertainty rows')
    result.loc[~eligible,'guide_ci_status']='not_in_original_plotted_cohort'
    result.to_csv(out/'guide_target_RNA_precision.csv',index=False)
    measured=result.loc[eligible&result.NTC_detected_fraction.ge(.1)]
    summary=dict(dataset=d['id'],cells=len(cells),NTC_cells=int(ntc.sum()),original_plotted_guides=int(eligible.sum()),
                 NTC_detection_ge10_guides=len(measured),estimated_CI=int(measured.guide_ci_status.eq('estimated').sum()),
                 all_zero_no_CI=int(measured.guide_ci_status.eq('all_zero_no_RNA_mean_interval').sum()),
                 original_points_unchanged=True,source_47GB_opened=False,new_public_cache_bytes=consumed_bytes,
                 resamples=config['resamples'],bootstrap_method='SciPy BCa; NTC point mean fixed',
                 NTC_range='same n without replacement; reference excludes selected NTC cells',
                 biological_replicates=False,significance_test=False,editing_efficiency=False,
                 python_versions=dict(numpy=np.__version__,scipy=scipy.__version__))
    atomic(out/'summary.json',summary);atomic(out/'input_proof.json',dict(contract=contract,source_contract=source_contract,**proof))
    atomic(out/'complete.json',dict(contract=contract,result=summary,files={n:sha(out/n) for n in ['guide_target_RNA_precision.csv','summary.json','input_proof.json']}))
    return summary


def run(config,resume=False):
    root=Path(config['output'])
    if root.parent!=Path(config['allowed_output_parent']):raise ValueError('Wrong output parent')
    root.mkdir(exist_ok=resume);lock=root/'owner.lock'
    owner=dict(pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19])
    with lock.open('x') as f:json.dump(owner,f)
    started=time.monotonic()
    def progress(stage,**kw):
        atomic(root/'status.json',dict(stage=stage,time=time.time(),elapsed=time.monotonic()-started,**owner,**kw));print(stage,kw,flush=True)
    def alarm(*_):raise TimeoutError('Precision time guard; target blocks and completed intervals preserved')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(config['max_seconds'])
    try:
        results={}
        for d in config['datasets']:
            out=root/d['id'];completed=out/'complete.json'
            if completed.exists():
                saved=json.loads(completed.read_text())
                if saved['contract']!=contracts(d,config)[1]:raise ValueError('Completed precision configuration differs')
                for name,h in saved['files'].items():
                    if sha(out/name)!=h:raise ValueError('Completed dataset changed')
                results[d['id']]=saved['result'];continue
            progress('prepare_intervals',dataset=d['id'])
            results[d['id']]=dataset(d,config,out,lambda s,**kw:progress(s,dataset=d['id'],**kw))
        atomic(root/'summary.json',dict(code_commit=config['expected_commit'],datasets=results))
        progress('complete_pending_combination_figure')
    except BaseException as exc:progress('failed_preserved',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    run(json.loads(a.config.read_text()),a.resume)
