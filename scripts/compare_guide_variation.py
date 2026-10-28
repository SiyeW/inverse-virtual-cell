"""Equal-cell-number RNA variation, retaining the original expression recipe.

Distances compare log2(1 + mean normalized RNA), not a correlation replot.
A common NTC reference cancels in pairwise differences. Sampled cell groups
are disjoint within a guide; repeated samples are NOT biological replicates.
"""
from __future__ import annotations
from pathlib import Path
import argparse,gzip,hashlib,io,json,os,signal,time
import numpy as np
import pandas as pd
from scipy import sparse
from compare_established_perturbseq import load_counts
from compare_agreement_overlay import members,fixed_plan,sha,atomic
from perturbseq_core_effects import stable_rng,pair_effects
from public_count_cache import check_arrays,axis_hash

CATEGORIES=['within_guide','same_target','different_target','NTC']


def bind_inputs(plan,cells,scale,key):
    contract=hashlib.sha256(plan['contract'].encode())
    for values in [cells.cell_barcode.astype(str),cells.guide_id.astype(str)]:contract.update(axis_hash(values).encode())
    contract.update(cells.total_counts.to_numpy(dtype='<f8').tobytes())
    contract.update(json.dumps(dict(scale=float(scale),source_key=key),sort_keys=True).encode())
    plan['contract']=contract.hexdigest()


def make_plan(cells,records,pairs,support,genes,common,n=10,repeats=8,seed=20261029):
    """Plan using identities and coverage only, before reading expression."""
    selected=records.loc[records.source_control_type.isin(['targeting','non_targeting_control'])&records.cells.ge(2*n)].sort_values('guide').reset_index(drop=True)
    eligible=set(selected.guide)
    pairs=pairs.loc[pairs.guide_a.isin(eligible)&pairs.guide_b.isin(eligible)].copy()
    plan=fixed_plan(selected,support,pairs,seed)
    if not len(plan):raise ValueError('No same-target pair with disjoint equal-size groups')
    guides=cells.guide_id.astype(str).to_numpy();positions={g:np.flatnonzero(guides==g) for g in selected.guide}
    for row in selected.itertuples():
        if len(positions[row.guide])!=int(row.cells):raise ValueError('Original guide cell count differs')
    sizes=[n//2,n];sizes=sorted(set(x for x in sizes if x>=2))
    rows=[];cols=[];weights=[];draws=[];row_lookup={};count=0
    for row in selected.itertuples():
        available=positions[row.guide]
        for repeat in range(repeats):
            chosen=stable_rng(seed,f'variation:{row.guide}:{repeat}').choice(available,2*n,replace=False)
            for size in sizes:
                for side in [0,1]:
                    use=chosen[side*n:side*n+size]
                    row_lookup[row.guide,size,repeat,side]=count
                    rows.extend([count]*size);cols.extend(use.tolist());weights.extend([1/size]*size)
                    draws.append(dict(guide=row.guide,cells=size,repeat=repeat+1,side=side,operator_row=count,analysis_rows=';'.join(map(str,use))))
                    count+=1
    full_start=count
    for row in selected.itertuples():
        use=positions[row.guide];rows.extend([count]*len(use));cols.extend(use.tolist());weights.extend([1/len(use)]*len(use));count+=1
    ntc=np.flatnonzero(cells.role.eq('non_targeting_control').to_numpy())
    if len(ntc)<2:raise ValueError('No NTC reference')
    reference_row=count;rows.extend([count]*len(ntc));cols.extend(ntc.tolist());weights.extend([1/len(ntc)]*len(ntc));count+=1
    operator=sparse.csr_matrix((weights,(rows,cols)),shape=(count,len(cells)))
    if not np.allclose(np.asarray(operator.sum(1)).ravel(),1,rtol=1e-12,atol=1e-12):raise ValueError('Incorrect averaging operator')
    # The same per-pair RNA mask is used for all four comparison categories.
    common_lookup={g:i for i,g in enumerate(common)};direct={}
    for row in records.itertuples():
        if pd.notna(row.direct_feature_index):
            g=str(genes[int(row.direct_feature_index)])
            if g in common_lookup:direct[row.guide]=common_lookup[g]
    keep=np.ones((len(plan),len(common)),bool)
    for i,row in enumerate(plan.itertuples()):
        for g in [row.guide_a,row.guide_b,row.other_a,row.other_b]:
            if g in direct:keep[i,direct[g]]=False
    ops=[]
    for size in sizes:
        comparisons=[]
        for left,right,sl,sr in [('guide_a','guide_a',0,1),('guide_b','guide_b',0,1),('guide_a','guide_b',0,0),('other_a','other_b',0,0),('NTC_a','NTC_b',0,0)]:
            a=np.array([[row_lookup[getattr(row,left),size,r,sl] for r in range(repeats)] for row in plan.itertuples()])
            b=np.array([[row_lookup[getattr(row,right),size,r,sr] for r in range(repeats)] for row in plan.itertuples()])
            comparisons.append((a,b))
        ops.append(comparisons)
    # Examples are chosen by coverage, not RNA differences/correlations.
    ordered=plan.assign(min_coverage=np.minimum(plan.cells_a,plan.cells_b)).sort_values(['min_coverage','target','guide_a','guide_b'],ascending=[False,True,True,True])
    case_i=int(ordered.index[0]);example=plan.iloc[case_i]
    case_ops=[(row_lookup[example.guide_a,n,0,0],row_lookup[example.guide_b,n,0,0]),
              (row_lookup[example.other_a,n,0,0],row_lookup[example.other_b,n,0,0])]
    contract=hashlib.sha256()
    for a in [operator.indptr,operator.indices,operator.data,keep]:contract.update(np.asarray(a).tobytes())
    contract.update(json.dumps(list(map(str,common))).encode())
    return dict(operator=operator,records=selected,plan=plan,draws=pd.DataFrame(draws),sizes=sizes,
        full_start=full_start,reference_row=reference_row,keep=keep,ops=ops,case_index=case_i,case_ops=case_ops,
        samples=repeats,ntc_cells=len(ntc),contract=contract.hexdigest(),n=n)


class Accumulator:
    def __init__(self,plan,expected,reference,checkpoint):
        self.plan=plan;self.expected=expected;self.reference=reference;self.checkpoint=Path(checkpoint)
        self.sums=np.zeros((len(plan['sizes']),5,len(plan['plan']),plan['samples']),float)
        self.done=np.zeros(expected.shape[1],bool);self.x=np.full((2,expected.shape[1]),np.nan);self.y=self.x.copy();self.max_error=0.
        receipt=self.checkpoint.with_suffix('.json')
        if receipt.exists():
            saved=json.loads(receipt.read_text())
            if saved['contract']!=plan['contract'] or sha(self.checkpoint)!=saved['sha256']:raise ValueError('Variation checkpoint contract/content differs')
            with np.load(self.checkpoint,allow_pickle=False) as z:
                self.sums=z['sums'];self.done=z['done'];self.x=z['x'];self.y=z['y'];self.max_error=float(z['max_error'])
            if self.sums.shape!=(len(plan['sizes']),5,len(plan['plan']),plan['samples']) or self.done.shape!=(expected.shape[1],):raise ValueError('Checkpoint shape differs')

    def consume(self,columns,means):
        columns=np.asarray(columns,int)
        if self.done[columns].any():raise ValueError('Completed RNA block must not be replayed')
        p=self.plan
        error=max(float(np.max(np.abs(means[p['full_start']:p['reference_row']]-self.expected[:,columns]))),
                  float(np.max(np.abs(means[p['reference_row']]-self.reference[columns]))))
        if error>1e-8:raise ValueError('Original full-guide mean/reference bridge failed: '+str(error))
        self.max_error=max(self.max_error,error)
        logged=np.log2(means+1);keep=p['keep'][:,columns]
        for s,comparisons in enumerate(p['ops']):
            for k,(a,b) in enumerate(comparisons):
                delta=np.abs(logged[a]-logged[b])
                self.sums[s,k]+=(delta*keep[:,None,:]).sum(axis=2)
        center=logged[p['reference_row']]
        for i,(a,b) in enumerate(p['case_ops']):
            self.x[i,columns]=logged[a]-center;self.y[i,columns]=logged[b]-center
        self.done[columns]=True
        tmp=self.checkpoint.with_suffix('.tmp')
        with tmp.open('wb') as f:np.savez_compressed(f,sums=self.sums,done=self.done,x=self.x,y=self.y,max_error=np.array(self.max_error))
        tmp.replace(self.checkpoint)
        atomic(self.checkpoint.with_suffix('.json'),dict(contract=p['contract'],sha256=sha(self.checkpoint),completed_RNA=int(self.done.sum()),mean_bridge_max_error=self.max_error))

    def finish(self,dataset,genes,names):
        if not self.done.all():raise ValueError('Incomplete variation RNA coverage')
        p=self.plan;count=p['keep'].sum(1)
        distance=self.sums/count[None,None,:,None]
        values=np.stack([(distance[:,0]+distance[:,1])/2,distance[:,2],distance[:,3],distance[:,4]],axis=1)
        rows=[];summaries=[]
        for s,size in enumerate(p['sizes']):
            for k,kind in enumerate(CATEGORIES):
                for i,q in enumerate(p['plan'].itertuples()):
                    for r in range(p['samples']):rows.append(dict(dataset=dataset,cells_per_group=size,comparison=kind,pair_id=i,target=q.target,repeat=r+1,mean_absolute_difference=float(values[s,k,i,r]),RNA_features=int(count[i])))
                    summaries.append(dict(dataset=dataset,cells_per_group=size,comparison=kind,pair_id=i,target=q.target,guide_a=q.guide_a,guide_b=q.guide_b,mean_absolute_difference=float(np.median(values[s,k,i])),RNA_features=int(count[i])))
        q=p['plan'].iloc[p['case_index']];case=[]
        for k,(kind,a,b) in enumerate([('same_target',q.guide_a,q.guide_b),('different_target',q.other_a,q.other_b)]):
            for j in np.flatnonzero(p['keep'][p['case_index']]):
                case.append(dict(dataset=dataset,comparison=kind,guide_a=a,guide_b=b,cells_per_group=p['n'],repeat=1,gene_id=str(genes[j]),gene_name=str(names[j]),x=float(self.x[k,j]),y=float(self.y[k,j])))
        return pd.DataFrame(rows),pd.DataFrame(summaries),pd.DataFrame(case)


def prepare_expected(baseline,records,common):
    saved,proof=members(baseline/'tables/rna_effect_profiles.npz',['means','reference_means','gene_ids','guides'])
    original=pd.read_csv(baseline/'tables/guide_effect_metadata.csv')
    if not np.array_equal(saved['guides'].astype(str),original.guide.astype(str)):raise ValueError('Original guide axes differ')
    columns=pd.Index(saved['gene_ids'].astype(str)).get_indexer(common)
    if (columns<0).any():raise ValueError('Common RNA not in original result')
    rows=pd.Index(saved['guides'].astype(str)).get_indexer(records.guide.astype(str))
    first=int(np.flatnonzero(original.source_control_type.ne('non_targeting_control'))[0])
    return saved['means'][rows][:,columns].astype(float),saved['reference_means'][first,columns].astype(float),saved['gene_ids'].astype(str),proof


def process_lab(cfg,common,out,n,repeats,seed,progress):
    baseline=Path(cfg['baseline']);records=pd.read_csv(baseline/'tables/guide_effect_metadata.csv')
    x,cells,genes,count_proof=load_counts(cfg['count_input'],progress)
    old,oldproof=members(baseline/'tables/rna_effect_profiles.npz',['gene_ids'])
    pairs=pd.read_csv(baseline/'tables/within_target_pairs.csv');support=pd.read_csv(baseline/'tables/target_support.csv')
    plan=make_plan(cells,records,pairs,support,old['gene_ids'].astype(str),common,n,repeats,seed)
    bind_inputs(plan,cells,cfg['scale'],count_proof.get('canonical_counts_sha256',count_proof.get('source_sha256')))
    expected,ref,_,proof=prepare_expected(baseline,plan['records'],common)
    ratio=float(cfg['scale'])/float(cfg.get('baseline_scale',cfg['scale']))
    expected*=ratio;ref*=ratio
    selected=pd.Index(genes.gene_id.astype(str)).get_indexer(common)
    if (selected<0).any():raise ValueError('Common RNA missing from count input')
    panel=x[:,selected].astype(np.float64).tocsc();del x
    factors=cfg['scale']/cells.total_counts.to_numpy(float)
    for first in range(0,len(panel.data),1024**2):
        last=min(len(panel.data),first+1024**2);panel.data[first:last]*=factors[panel.indices[first:last]]
    accumulator=Accumulator(plan,expected,ref,out/'rna_checkpoint.npz')
    for first in range(0,len(common),64):
        columns=np.arange(first,min(first+64,len(common)))
        if accumulator.done[columns].all():continue
        if accumulator.done[columns].any():raise ValueError('Partial lab RNA span cannot be silently changed')
        means=(plan['operator']@panel[:,columns]).toarray();accumulator.consume(columns,means)
        progress('variation_RNA_blocks',completed_RNA=int(accumulator.done.sum()),total_RNA=len(common))
    names=genes.gene_name.to_numpy(str)[selected]
    return plan,accumulator,cells,names,dict(counts=count_proof,baseline_members=proof,mean_bridge_max_error=accumulator.max_error,scale=cfg['scale'],original_scale=cfg.get('baseline_scale',cfg['scale']),linear_bridge_scaling=ratio)


def process_public(cfg,common,out,n,repeats,seed,progress):
    baseline=Path(cfg['baseline']);manifest=json.loads((baseline/'raw_count_cache/manifest.json').read_text())
    if sha(baseline/'raw_count_cache/manifest.json')!=cfg['cache_manifest_sha256']:raise ValueError('Public count cache manifest changed')
    all_obs=pd.read_csv(baseline/'tables/all_released_cells_with_additional_QC.csv.gz',index_col=0)
    genes=pd.read_csv(baseline/'tables/all_RNA_feature_selection.csv.gz',index_col=0)
    if axis_hash(all_obs.index)!=manifest['contract']['rows'] or axis_hash(genes.index)!=manifest['contract']['genes']:raise ValueError('Public cache axes differ')
    cells=all_obs.loc[all_obs.main_analysis_stage1.astype(bool)].copy();cells['cell_barcode']=cells.index.astype(str)
    cells['role']=np.where(cells.gene_target.eq('NTC'),'non_targeting_control','targeting')
    records=pd.read_csv(baseline/'tables/guide_effect_metadata.csv');pairs=pd.read_csv(baseline/'tables/within_target_pairs.csv');support=pd.read_csv(baseline/'tables/target_support.csv')
    original,_=members(baseline/'tables/rna_effect_profiles.npz',['gene_ids'])
    plan=make_plan(cells,records,pairs,support,original['gene_ids'].astype(str),common,n,repeats,seed)
    bind_inputs(plan,cells,cfg['scale'],cfg['cache_manifest_sha256'])
    expected,ref,_,proof=prepare_expected(baseline,plan['records'],common)
    ratio=float(cfg['scale'])/float(cfg.get('baseline_scale',cfg['scale']))
    expected*=ratio;ref*=ratio
    accumulator=Accumulator(plan,expected,ref,out/'rna_checkpoint.npz')
    columns=genes.index.astype(str).get_indexer(common)
    if (columns<0).any():raise ValueError('Public common RNA missing')
    remap=np.full(len(all_obs),-1,np.int32);remap[cells.source_row.to_numpy(int)]=np.arange(len(cells),dtype=np.int32)
    factors=cfg['scale']/cells.total_counts.to_numpy(float);consumed=[]
    for entry in manifest['blocks']:
        use=np.flatnonzero((columns>=entry['first'])&(columns<entry['end']))
        if not len(use) or accumulator.done[use].all():continue
        if accumulator.done[use].any():raise ValueError('Partially accepted public source block')
        path=baseline/'raw_count_cache'/entry['file'];before=path.stat()
        if before.st_size!=entry['bytes']:raise ValueError('Public block size changed')
        payload=path.read_bytes()
        if len(payload)!=entry['bytes'] or hashlib.sha256(payload).hexdigest()!=entry['sha256']:raise ValueError('Public accepted count block changed')
        with np.load(io.BytesIO(payload),allow_pickle=False) as z:data,rows,ptr=z['data'],z['rows'],z['ptr']
        check_arrays(data,rows,ptr,len(all_obs));values=[];indices=[];pointers=[0]
        for position in use:
            j=columns[position]-entry['first'];a,b=int(ptr[j]),int(ptr[j+1]);r=remap[rows[a:b]];keep=r>=0;rr=r[keep]
            values.append(data[a:b][keep].astype(float)*factors[rr]);indices.append(rr);pointers.append(pointers[-1]+len(rr))
        matrix=sparse.csc_matrix((np.concatenate(values),np.concatenate(indices),pointers),shape=(len(cells),len(use)))
        if (path.stat().st_size,path.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):raise ValueError('Public block changed while reading')
        means=(plan['operator']@matrix).toarray();accumulator.consume(use,means)
        consumed.append(dict(file=entry['file'],bytes=entry['bytes'],sha256=entry['sha256']))
        progress('variation_public_cache',completed_RNA=int(accumulator.done.sum()),total_RNA=len(common),cached_bytes=sum(v['bytes'] for v in consumed))
    names=genes.gene_name.to_numpy(str)[columns]
    return plan,accumulator,cells,names,dict(source_count_payload_read=False,accepted_cache_blocks_read=consumed,baseline_members=proof,mean_bridge_max_error=accumulator.max_error,scale=cfg['scale'],original_scale=cfg.get('baseline_scale',cfg['scale']),linear_bridge_scaling=ratio)


def run(config,resume=False):
    root=Path(config['output'])
    if not root.is_absolute() or root.parent!=Path(config['allowed_output_parent']):raise ValueError('Explicit independent output parent required')
    root.mkdir(parents=True,exist_ok=resume)
    if (root/'owner.lock').exists():raise ValueError('Existing variation owner')
    owner=dict(pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19])
    lock=root/'owner.lock'
    with lock.open('x') as stream:json.dump(owner,stream)
    started=time.monotonic()
    def progress(stage,**kw):
        atomic(root/'status.json',dict(stage=stage,**owner,elapsed=time.monotonic()-started,time=time.time(),**kw));print(stage,kw,flush=True)
    def alarm(*_):raise TimeoutError('Finite variation protection reached; dataset/RNA checkpoints preserved')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(config.get('max_seconds',5400))
    try:
        common=pd.read_csv(config['common_RNA']).gene_id.astype(str).to_numpy()
        if sha(Path(config['common_RNA']))!=config['common_RNA_sha256'] or len(common)!=2272:raise ValueError('Fixed common RNA changed')
        results={}
        for d in config['datasets']:
            name=d['id'];out=root/name;out.mkdir(exist_ok=resume)
            for filename,expected_sha in d.get('metadata_sha256',{}).items():
                if sha(Path(filename))!=expected_sha:raise ValueError('Frozen metadata changed: '+filename)
            complete=out/'complete.json'
            if complete.exists():
                entry=json.loads(complete.read_text())
                for filename,row in entry['files'].items():
                    if sha(out/filename)!=row['sha256']:raise ValueError('Completed variation artifact changed')
                results[name]=entry['result'];continue
            def update(stage,**kw):progress(stage,dataset=name,**kw)
            update('prepare_variation')
            function=process_public if d['kind']=='public_cache' else process_lab
            plan,acc,cells,names,proof=function(d,common,out,config['n'],config['repeats'],config['seed'],update)
            plan['records'].to_csv(out/'eligible_guides.csv',index=False);plan['plan'].to_csv(out/'comparison_plan.csv',index=False)
            plan['draws'].to_csv(out/'sampling_plan.csv.gz',index=False)
            draws,summary,cases=acc.finish(name,common,names)
            draws.to_csv(out/'distance_draws.csv.gz',index=False);summary.to_csv(out/'variation_points.csv',index=False);cases.to_csv(out/'RNA_case_points.csv.gz',index=False)
            row=plan['plan'].iloc[plan['case_index']]
            metadata=dict(query_target=row.target,guide_a=row.guide_a,guide_b=row.guide_b,other_a=row.other_a,other_b=row.other_b,
                          full_cells_a=int(row.cells_a),full_cells_b=int(row.cells_b),cells_per_group=config['n'],repeat=1,
                          selection='largest minimum coverage, then lexical target/guide; no expression outcome selection')
            atomic(out/'case_selection.json',metadata);atomic(out/'input_proof.json',proof)
            result=dict(analysis_cells=len(cells),NTC_cells=plan['ntc_cells'],eligible_guides=len(plan['records']),same_target_pairs=len(plan['plan']),targets=int(plan['plan'].target.nunique()),n=config['n'],repeats=config['repeats'],RNA_base=len(common),bridge_max_error=acc.max_error,case=metadata)
            files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file() and p.name not in ['complete.json','rna_checkpoint.npz','rna_checkpoint.json']}
            atomic(complete,dict(files=files,result=result));results[name]=result;update('dataset_complete')
        summary=dict(complete=True,datasets=results,code_commit=config['expected_commit'],method='per-cell total RNA normalization to 10000; log2(mean+1); same reference cancels in pair distance',
                     distance='mean absolute RNA difference',sizes=[config['n']//2,config['n']],primary_n=config['n'],repeats=config['repeats'],
                     within_guide='average of the two query guides disjoint split distances',RNA_mask='union of direct target RNAs for the query and different-target comparator, applied to all four arms',
                     batch_conditioning=False,cell_filtering_changed=False,biological_replicates=False,editing_efficiency_assessed=False)
        atomic(root/'summary.json',summary);progress('complete_pending_plot_review')
    except BaseException as exc:progress('failed_preserved',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--resume',action='store_true');args=p.parse_args()
    run(json.loads(args.config.read_text()),args.resume)
