"""Fixed-cohort old/new agreement comparison for lab and public data.

The common revision changes only mean/log order. References remain pooled NTC
in BOTH methods and ALL datasets; this does not pretend to batch-correct an
unlabelled laboratory matrix. No cell filtering, new gene selection or hit test.
"""
from __future__ import annotations
import argparse,gc,hashlib,io,json,os,signal,time,zipfile
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from compare_established_perturbseq import load_counts
from perturbseq_core_effects import pair_effects,stable_rng


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()


def atomic(p,o):
    q=p.with_suffix('.tmp');q.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n');q.replace(p)


def members(path,names):
    """Read each requested NPZ member once, CRC checked, not the whole archive."""
    out={};receipt=[]
    before=Path(path).stat()
    with zipfile.ZipFile(path) as z:
        for name in names:
            entry=z.getinfo(name+'.npy')
            if entry.file_size>100*1024**2:raise ValueError('Unexpected profile member size')
            body=z.read(entry.filename)
            out[name]=np.load(io.BytesIO(body),allow_pickle=False)
            receipt.append(dict(member=entry.filename,compressed_bytes=entry.compress_size,
                                bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),crc_verified=True))
    after=Path(path).stat()
    if (before.st_size,before.st_mtime_ns,before.st_dev)!=(after.st_size,after.st_mtime_ns,after.st_dev):raise ValueError('Profile identity changed')
    return out,dict(path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,consumed_members=receipt,whole_file_rehashed=False)


def check_axes(records,guides,genes):
    if not records.guide.is_unique or not pd.Index(genes).is_unique:raise ValueError('Duplicate profile identity')
    if not np.array_equal(records.guide.astype(str),guides.astype(str)):raise ValueError('Guide profile order differs')


def centers_from_means(means,records,ntc_mean,ntc_n):
    refs=np.broadcast_to(ntc_mean,means.shape).copy()
    nc=records.source_control_type.eq('non_targeting_control').to_numpy()
    for i in np.flatnonzero(nc):
        n=int(records.cells.iloc[i])
        if ntc_n-n<2:raise ValueError('Insufficient leave-whole-guide-out NTC reference')
        refs[i]=(ntc_n*ntc_mean-n*means[i])/(ntc_n-n)
    return refs


def apply_mask(values,records):
    values=values.copy()
    for i,row in enumerate(records.itertuples()):
        if pd.notna(row.direct_feature_index):values[i,int(row.direct_feature_index)]=np.nan
    return values


def linear_and_log_means(x,cells,records,scale):
    """Same selected cells/genes, double precision, whole-guide NTC exclusion."""
    labels=cells.guide_id.astype(str).to_numpy();lookup={g:i for i,g in enumerate(records.guide.astype(str))}
    nc=cells.role.eq('non_targeting_control').to_numpy();n=int(nc.sum())
    if n<2:raise ValueError('No NTC cells')
    totals=cells.total_counts.to_numpy(float)
    if (totals<=0).any():raise ValueError('Nonpositive RNA total')
    row_index=[];col_index=[];weights=[]
    counts=pd.Series(labels).value_counts()
    for r in records.itertuples():
        if int(counts.get(r.guide,0))!=int(r.cells):raise ValueError('Original guide cell count changed')
    for i,g in enumerate(labels):
        if g in lookup:
            row_index.append(lookup[g]);col_index.append(i);weights.append(1/float(counts[g]))
        if nc[i]:row_index.append(len(records));col_index.append(i);weights.append(1/n)
    operator=sparse.csr_matrix((weights,(row_index,col_index)),shape=(len(records)+1,len(cells)))
    norm=x.astype(np.float64)
    for first in range(0,len(cells),256):
        last=min(first+256,len(cells));a,b=norm.indptr[first],norm.indptr[last]
        norm.data[a:b]*=np.repeat(scale/totals[first:last],np.diff(norm.indptr[first:last+1]))
    linear=(operator@norm).toarray()
    np.add(norm.data,1,out=norm.data);np.log2(norm.data,out=norm.data)
    logged=(operator@norm).toarray()
    del norm
    return linear[:-1],linear[-1],logged[:-1],logged[-1],n


def lab_profiles(config,baseline,progress):
    records=pd.read_csv(baseline/'tables/guide_effect_metadata.csv')
    old,proof=members(baseline/'tables/rna_effect_profiles.npz',['means','reference_means','gene_ids','guides'])
    check_axes(records,old['guides'],old['gene_ids'])
    x,cells,genes,count_proof=load_counts(config['count_input'],progress)
    ids=pd.Index(genes.gene_id.astype(str));selected=ids.get_indexer(old['gene_ids'].astype(str))
    if (selected<0).any():raise ValueError('Original selected RNA absent from counts')
    # Keep original per-cell totals for scaling; do not divide by selected genes.
    chosen=x[:,selected].tocsr();del x;gc.collect()
    lm,lc,bm,bc,n=linear_and_log_means(chosen,cells,records,float(config['scale']))
    del chosen;gc.collect()
    ar=centers_from_means(lm,records,lc,n)
    mean_error=float(np.max(np.abs(lm-old['means'])))
    ref_error=float(np.max(np.abs(ar-old['reference_means'])))
    if max(mean_error,ref_error)>1e-8:raise ValueError('Original all-cell mean/reference bridge failed: '+str((mean_error,ref_error)))
    old_effect=apply_mask(np.log2((lm+1)/(ar+1)),records)
    new_effect=apply_mask(bm-centers_from_means(bm,records,bc,n),records)
    return records,old['gene_ids'].astype(str),old_effect,new_effect,dict(
        counts=count_proof,baseline_members=proof,mean_bridge_max_error=mean_error,
        reference_bridge_max_error=ref_error,NTC_cells=n,analysis_cells=len(cells),normalization_target_umi=config['scale'])


def public_profiles(config,baseline,progress):
    records=pd.read_csv(baseline/'tables/guide_effect_metadata.csv')
    old,proof=members(baseline/'tables/rna_effect_profiles.npz',['means','reference_means','gene_ids','guides'])
    check_axes(records,old['guides'],old['gene_ids'])
    revision=Path(config['revision']);meta=pd.read_csv(revision/'tables/group_metadata.csv')
    logged,logproof=members(revision/'tables/B1_within_batch_cellwise_log2_profiles.npz',['means','keys','gene_ids'])
    if not np.array_equal(logged['keys'].astype(str),meta.key.astype(str)) or not np.array_equal(logged['gene_ids'].astype(str),old['gene_ids'].astype(str)):raise ValueError('Revision axes differ')
    bm=np.empty_like(old['means']);n_total=0;ntc_sum=np.zeros(len(old['gene_ids']))
    for i,r in enumerate(records.itertuples()):
        selected=meta.guide.eq(r.guide).to_numpy();w=meta.loc[selected,'cells'].to_numpy(float)
        if int(w.sum())!=int(r.cells):raise ValueError('Public guide count changed')
        bm[i]=np.average(logged['means'][selected],axis=0,weights=w)
        if r.source_control_type=='non_targeting_control':n_total+=int(r.cells);ntc_sum+=int(r.cells)*bm[i]
    if n_total!=config['NTC_cells']:raise ValueError('Public pooled controls differ')
    nc_mean=ntc_sum/n_total
    old_effect=apply_mask(np.log2((old['means']+1)/(old['reference_means']+1)),records)
    new_effect=apply_mask(bm-centers_from_means(bm,records,nc_mean,n_total),records)
    return records,old['gene_ids'].astype(str),old_effect,new_effect,dict(
        baseline_members=proof,revision_members=logproof,source_count_payload_read=False,
        NTC_cells=n_total,analysis_cells=int(records.cells.sum()),normalization_target_umi=config['scale'],
        note='Only saved per-cell-log means reused. Within-batch references/effects NOT reused in this common comparison.')


def fixed_plan(records,support,pairs,seed=20261029):
    """Freeze one ordinary-target cohort and controls for both methods."""
    target_roles=support.set_index('target').source_control_type
    ordinary=pairs.loc[pairs.target.map(target_roles).eq('targeting')].copy()
    ordinary=ordinary.loc[ordinary.correlation_defined.astype(bool)]
    lookup=records.set_index('guide');tc=records.loc[records.source_control_type.eq('targeting')]
    nc=records.loc[records.source_control_type.eq('non_targeting_control')]
    ncpairs=list(combinations(nc.guide.astype(str),2))
    if not ncpairs:raise ValueError('No control guide pairs')
    nclog=np.log(np.array([[lookup.loc[a,'cells'],lookup.loc[b,'cells']] for a,b in ncpairs],float))
    rows=[];guides=tc.guide.astype(str).to_numpy();targets=tc.target.astype(str).to_numpy();sizes=tc.cells.to_numpy(float)
    for row in ordinary.itertuples():
        target_n=np.log([lookup.loc[row.guide_a,'cells'],lookup.loc[row.guide_b,'cells']])
        rng=stable_rng(seed,row.guide_a+'|'+row.guide_b)
        ia=rng.integers(0,len(tc),4096);ib=rng.integers(0,len(tc),4096)
        keep=(targets[ia]!=targets[ib])&(targets[ia]!=str(row.target))&(targets[ib]!=str(row.target))
        ia,ib=ia[keep],ib[keep]
        if not len(ia):raise ValueError('No ordinary different-target control')
        score=(np.log(sizes[ia])-target_n[0])**2+(np.log(sizes[ib])-target_n[1])**2
        j=int(np.argmin(score));other=(guides[ia[j]],guides[ib[j]])
        forward=((nclog-target_n)**2).sum(1);backward=((nclog[:,::-1]-target_n)**2).sum(1)
        j=int(np.argmin(np.minimum(forward,backward)));control=ncpairs[j]
        if backward[j]<forward[j]:control=control[::-1]
        rows.append(dict(target=row.target,guide_a=row.guide_a,guide_b=row.guide_b,
            other_a=other[0],other_b=other[1],NTC_a=control[0],NTC_b=control[1],
            cells_a=int(np.exp(target_n[0])+.5),cells_b=int(np.exp(target_n[1])+.5),
            other_cells_a=int(lookup.loc[other[0],'cells']),other_cells_b=int(lookup.loc[other[1],'cells']),
            NTC_cells_a=int(lookup.loc[control[0],'cells']),NTC_cells_b=int(lookup.loc[control[1],'cells'])))
    return pd.DataFrame(rows)


def evaluate(records,effects,plan,support,method):
    lookup={g:i for i,g in enumerate(records.guide.astype(str))};rows=[]
    for row in plan.itertuples():
        for kind,a,b in [('same_target',row.guide_a,row.guide_b),('different_target',row.other_a,row.other_b),('NTC',row.NTC_a,row.NTC_b)]:
            value=pair_effects(effects[lookup[a]],effects[lookup[b]])
            if not value['correlation_defined']:raise ValueError('Undefined fixed comparison')
            rows.append(dict(method=method,target=row.target,comparison=kind,guide_a=a,guide_b=b,
                             pearson=value['pearson'],features=value['features']))
    pair=pd.DataFrame(rows)
    med=pair.loc[pair.comparison.eq('same_target')].groupby('target').pearson.median()
    wanted=support.loc[support.source_control_type.eq('targeting')&support.pairs.gt(0),['target','largest_guide_fraction']].copy()
    if set(wanted.target)!=set(med.index):raise ValueError('Fixed target-point coverage changed')
    wanted['median_pair_r']=wanted.target.map(med);wanted['method']=method
    return wanted,pair


def run(config):
    out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
    lock=out/'owner.lock';lock.write_text(json.dumps(dict(pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19])))
    start=time.monotonic()
    def progress(stage,**kw):
        atomic(out/'status.json',dict(stage=stage,pid=os.getpid(),time=time.time(),elapsed=time.monotonic()-start,**kw));print(stage,kw,flush=True)
    def timeout(*_):raise TimeoutError('Finite overlay run protection reached; completed profile checkpoints preserved')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(config.get('max_seconds',5400))
    try:
        all_data={};proofs={}
        for d in config['datasets']:
            name=d['id'];baseline=Path(d['baseline']);progress('fixed_profiles',dataset=name)
            function=public_profiles if d['kind']=='saved_public_means' else lab_profiles
            rec,genes,a,b,proof=function(d,baseline,lambda stage,**kw:progress(stage,dataset=name,**kw))
            np.savez_compressed(out/(name+'_profiles_checkpoint.npz'),old=a,new=b,gene_ids=np.asarray(genes,dtype=str),guides=rec.guide.to_numpy(dtype=str))
            rec.to_csv(out/(name+'_records.csv'),index=False);atomic(out/(name+'_input_proof.json'),proof)
            support=pd.read_csv(baseline/'tables/target_support.csv');pairs=pd.read_csv(baseline/'tables/within_target_pairs.csv')
            plan=fixed_plan(rec,support,pairs);plan.to_csv(out/(name+'_comparison_plan.csv'),index=False)
            old_points,old_pairs=evaluate(rec,a,plan,support,'old')
            joined=old_points.merge(support[['target','median_pair_r']],on='target',suffixes=('_recomputed','_baseline'))
            error=float(np.max(np.abs(joined.median_pair_r_recomputed-joined.median_pair_r_baseline)))
            if error>1e-8:raise ValueError('Original scatter parity failed '+str(error))
            proof['original_scatter_max_error']=error;atomic(out/(name+'_input_proof.json'),proof)
            all_data[name]=(rec,genes,a,b,support,plan)
            for method,effect in [('old',a),('new',b)]:
                points,pair=evaluate(rec,effect,plan,support,method)
                points.to_csv(out/f'{name}_{method}_points.csv',index=False);pair.to_csv(out/f'{name}_{method}_pairs.csv',index=False)
            proofs[name]=proof;progress('dataset_complete',dataset=name,targets=len(old_points))
            gc.collect()
        common=sorted(set.intersection(*(set(v[1]) for v in all_data.values())))
        if len(common)<100:raise ValueError('Insufficient common RNA for sensitivity')
        pd.DataFrame({'gene_id':common}).to_csv(out/'common_RNA.csv',index=False)
        for name,(rec,genes,a,b,support,plan) in all_data.items():
            index=pd.Index(genes).get_indexer(common)
            for method,effect in [('old',a[:,index]),('new',b[:,index])]:
                points,pair=evaluate(rec,effect,plan,support,method)
                points.to_csv(out/f'{name}_{method}_commonRNA_points.csv',index=False)
                pair.to_csv(out/f'{name}_{method}_commonRNA_pairs.csv',index=False)
        summary=dict(complete=True,method_old='log2((mean(normalized RNA)+1)/(mean(NTC)+1))',
                     method_new='mean(log2(1+normalized RNA))-mean_NTC(log2(1+normalized RNA))',
                     reference_both_methods='pooled NTC, whole own NTC guide excluded',batch_correction_in_main_comparison=False,
                     scientific_input_modified=False,new_cell_filtering=False,new_RNA_selection=False,
                     common_RNA_count=len(common),code_commit=config['expected_commit'],input_proofs=proofs,
                     inferential_pvalues=False,functional_hits=False,NTC_background='Nearest available guide-cell counts; residual size differences retained in plans. No false-positive-rate claim.')
        atomic(out/'summary.json',summary)
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file() and p.name not in ['owner.lock','status.json','output_manifest.json']}
        atomic(out/'output_manifest.json',manifest);progress('complete_pending_figure_review')
    except BaseException as exc:progress('failed_preserved',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    run(json.loads(args.config.read_text()))
