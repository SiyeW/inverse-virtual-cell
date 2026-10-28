"""Direct target RNA abundance; this does not estimate DNA editing efficiency.

Uses exact catalogue gene IDs, the frozen single-guide QC cohorts, all cells
per guide, full-cell library sizes, and pooled NTC for the shared comparison.
Public same-batch NTC values are supplementary, never substituted silently.
"""
from pathlib import Path
import argparse, hashlib, io, json, os, signal, time
import numpy as np
import pandas as pd
from scipy import sparse
from compare_established_perturbseq import load_counts
from compare_agreement_overlay import sha, atomic
from public_count_cache import axis_hash, check_arrays


def resolve_targets(catalogue, genes):
    """Do not infer host transcripts by parsing PB/ORF labels or gene symbols."""
    c=catalogue.copy()
    if not c.guide.is_unique or not pd.Index(genes.gene_id).is_unique:
        raise ValueError('Ambiguous guide or expression RNA ID')
    ids=pd.Index(genes.gene_id.astype(str))
    targeted=c.source_control_type.isin(['targeting','positive_control'])
    c['rna_column']=ids.get_indexer(c.target_gene_id.astype(str))
    c['mapping_status']=np.where(~targeted,'not_targeted',np.where(c.rna_column.ge(0),'exact_catalogue_id','missing_from_full_reference'))
    c.loc[~targeted,'rna_column']=-1
    return c


def summarize_gene(values, cells, design, gene_id, gene_name):
    """Zeros remain in each mean and detection denominator; no imputation."""
    values=np.asarray(values,float)
    if values.shape!=(len(cells),) or not np.isfinite(values).all() or (values<0).any():
        raise ValueError('Invalid normalized RNA values')
    ntc=cells.role.eq('non_targeting_control').to_numpy()
    if ntc.sum()<2:raise ValueError('Missing NTC reference')
    ref=float(values[ntc].mean()); detection=float((values[ntc]>0).mean())
    out=[]
    for row in design.itertuples():
        use=cells.guide_id.eq(row.guide).to_numpy();n=int(use.sum())
        v=values[use];mean=float(v.mean()) if n else np.nan
        ratio=mean/ref if n and ref>0 else np.nan
        r=dict(guide=row.guide,target_gene_id=gene_id,target=row.target_gene_name,
               measured_RNA_name=gene_name,source_control_type=row.source_control_type,
               cells=n,NTC_cells=int(ntc.sum()),mean_per_10000=mean,
               NTC_mean_per_10000=ref,detected_cells=int((v>0).sum()),
               detected_fraction=float((v>0).mean()) if n else np.nan,
               NTC_detected_fraction=detection,NTC_detected_cells=int((values[ntc]>0).sum()),
               RNA_remaining_ratio=ratio,reference_status='not_detected' if ref==0 else 'detected',
               mapping_status='exact_catalogue_id',zero_count_is_not_editing_proof=True)
        if 'batch' in cells:
            batches=cells.loc[use,'batch'].astype(str).unique()
            if n and len(batches)!=1:raise ValueError('Targeting guide spans multiple public batches')
            r['guide_batch']=batches[0] if n else ''
            refmask=ntc & cells.batch.astype(str).eq(batches[0]).to_numpy() if n else np.zeros(len(cells),bool)
            bm=float(values[refmask].mean()) if refmask.any() else np.nan
            r.update(same_batch_NTC_cells=int(refmask.sum()),same_batch_NTC_mean_per_10000=bm,
                     same_batch_RNA_remaining_ratio=mean/bm if n and bm>0 else np.nan)
        out.append(r)
    return pd.DataFrame(out)


def validate_cache(path, receipt, contract):
    r=json.loads(receipt.read_text())
    if r['contract']!=contract or sha(path)!=r['sha256']:
        raise ValueError('Target-count checkpoint identity/content differs')


def public_inputs(d):
    b=Path(d['baseline']);manifest=json.loads((b/'raw_count_cache/manifest.json').read_text())
    if sha(b/'raw_count_cache/manifest.json')!=d['cache_manifest_sha256']:raise ValueError('Public manifest changed')
    obs=pd.read_csv(b/'tables/all_released_cells_with_additional_QC.csv.gz',index_col=0)
    g=pd.read_csv(b/'tables/all_RNA_feature_selection.csv.gz',index_col=0)
    if axis_hash(obs.index)!=manifest['contract']['rows'] or axis_hash(g.index)!=manifest['contract']['genes']:
        raise ValueError('Public full RNA/cell axes changed')
    cells=obs.loc[obs.main_analysis_stage1.astype(bool)].copy()
    cells['cell_barcode']=cells.index.astype(str)
    cells['role']=np.where(cells.gene_target.eq('NTC'),'non_targeting_control','targeting')
    c=pd.read_csv(d['target_catalogue'],keep_default_na=False)
    c['guide']=c['id'].map(lambda s:'_'.join(s.rsplit('-',1)))
    if not c.guide.is_unique:raise ValueError('Ambiguous author guide conversion')
    c=c.loc[c.guide.isin(obs.guide_id)].copy()
    c['source_control_type']=np.where(c.target_gene_name.eq('Non-Targeting'),'non_targeting_control','targeting')
    lookup=c.set_index('guide').target_gene_name
    if not np.array_equal(cells.guide_id.map(lookup),cells.gene_target.replace({'NTC':'Non-Targeting'})):
        raise ValueError('Author target identities differ from cells')
    genes=g.reset_index().rename(columns={g.index.name or 'index':'gene_id'})
    return cells,genes,c,manifest,obs


def run_dataset(d,out,progress):
    out.mkdir(exist_ok=True)
    for p,h in d['metadata_sha256'].items():
        if sha(Path(p))!=h:raise ValueError('Frozen metadata changed: '+p)
    if sha(Path(d['target_catalogue']))!=d['target_catalogue_sha256']:raise ValueError('Target design changed')
    checkpoint=out/'target_counts.npz';receipt=out/'target_counts.json'
    contract=hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest()
    proof={}
    if d['kind']=='public_cache':
        cells,genes,c,manifest,all_obs=public_inputs(d)
    else:
        if receipt.exists():
            validate_cache(checkpoint,receipt,contract)
            cells=pd.read_csv(out/'cells.csv.gz',keep_default_na=False);genes=pd.read_csv(out/'genes.csv.gz',keep_default_na=False)
            for name,h in json.loads(receipt.read_text())['metadata'].items():
                if sha(out/name)!=h:raise ValueError('Cached metadata changed')
            proof=json.loads((out/'count_input_proof.json').read_text())
        else:
            x,cells,genes,proof=load_counts(d['count_input'],progress)
        c=pd.read_csv(d['target_catalogue'],keep_default_na=False)
        c=c[['guide','target_gene_id','target_gene_name','source_control_type']].copy()
    mapping=resolve_targets(c,genes);mapping.to_csv(out/'target_mapping.csv',index=False)
    valid=mapping.loc[mapping.rna_column.ge(0)]
    positions=np.sort(valid.rna_column.unique())
    if d['kind']!='public_cache':
        if receipt.exists():
            panel=sparse.load_npz(checkpoint)
        else:
            # Preserve the compact exact integer count panel before evaluation.
            panel=x[:,positions].tocsc();del x
            cells.to_csv(out/'cells.csv.gz',index=False);genes.to_csv(out/'genes.csv.gz',index=False)
            sparse.save_npz(checkpoint,panel)
            atomic(out/'count_input_proof.json',proof)
            atomic(receipt,dict(contract=contract,sha256=sha(checkpoint),positions=positions.tolist(),
                               metadata={n:sha(out/n) for n in ['cells.csv.gz','genes.csv.gz','count_input_proof.json']}))
        if panel.shape!=(len(cells),len(positions)) or json.loads(receipt.read_text())['positions']!=positions.tolist():
            raise ValueError('Cached target counts axes differ')
    cells['guide_id']=cells.guide_id.astype(str)
    if not cells.cell_barcode.is_unique or (cells.total_counts<=0).any():raise ValueError('Invalid cells/totals')
    factors=10000./cells.total_counts.to_numpy(float)
    chunks=[];consumed=[]
    if d['kind']=='public_cache':
        remap=np.full(len(all_obs),-1,np.int32);remap[cells.source_row.to_numpy(int)]=np.arange(len(cells))
        entries=manifest['blocks']
    else:entries=[dict(first=0,end=len(genes),file='lab_target_panel',bytes=checkpoint.stat().st_size)]
    for entry in entries:
        use=positions[(positions>=entry['first'])&(positions<entry['end'])]
        if not len(use):continue
        result=out/('block_%05d.csv'%entry['first']);rfile=result.with_suffix('.json')
        if rfile.exists():
            validate_cache(result,rfile,contract);chunks.append(pd.read_csv(result));continue
        if d['kind']=='public_cache':
            source=Path(d['baseline'])/'raw_count_cache'/entry['file'];before=source.stat()
            payload=source.read_bytes()
            if len(payload)!=entry['bytes'] or hashlib.sha256(payload).hexdigest()!=entry['sha256']:
                raise ValueError('Verified public count block changed')
            with np.load(io.BytesIO(payload),allow_pickle=False) as z:data,rows,ptr=z['data'],z['rows'],z['ptr']
            check_arrays(data,rows,ptr,len(all_obs));del payload
            if (source.stat().st_size,source.stat().st_mtime_ns)!=(before.st_size,before.st_mtime_ns):raise ValueError('Source changed')
            consumed.append(dict(file=entry['file'],bytes=entry['bytes'],sha256=entry['sha256']))
        frames=[]
        for j in use:
            if d['kind']=='public_cache':
                k=j-entry['first'];a,b=int(ptr[k]),int(ptr[k+1]);rr=remap[rows[a:b]];keep=rr>=0
                values=np.zeros(len(cells));values[rr[keep]]=data[a:b][keep]*factors[rr[keep]]
            else:
                k=int(np.searchsorted(positions,j));values=panel[:,k].toarray().ravel()*factors
            design=valid.loc[valid.rna_column.eq(j)]
            frames.append(summarize_gene(values,cells,design,str(genes.iloc[j].gene_id),str(genes.iloc[j].gene_name)))
        frame=pd.concat(frames,ignore_index=True);frame.to_csv(result,index=False)
        atomic(rfile,dict(contract=contract,sha256=sha(result),RNA_columns=use.tolist(),source=entry))
        chunks.append(frame);progress('target_RNA_blocks',completed_target_RNA=sum(x.target_gene_id.nunique() for x in chunks),total_target_RNA=len(positions))
    table=pd.concat(chunks,ignore_index=True)
    missing=mapping.loc[mapping.mapping_status.eq('missing_from_full_reference')].copy()
    if len(missing):
        missing=missing.rename(columns={'target_gene_name':'target'});missing['cells']=missing.guide.map(cells.guide_id.value_counts()).fillna(0).astype(int)
        table=pd.concat([table,missing],ignore_index=True)
    table.insert(0,'dataset',d['id']);table.to_csv(out/'guide_target_RNA.csv',index=False)
    atomic(out/'input_proof.json',dict(counts=proof,cache_blocks_read=consumed,source_47GB_opened=False,
         normalization='10000 * target integer count / full-cell RNA integer total',
         reference='all QC single-guide NTC cells, pooled for every primary dataset',
         targets='exact original catalogue gene_id; all RNA features, no response-gene filter',cells=len(cells),NTC_cells=int(cells.role.eq('non_targeting_control').sum())))
    eligible=table.loc[table.cells.ge(10)&table.mapping_status.eq('exact_catalogue_id')&table.NTC_mean_per_10000.gt(0)]
    measured=eligible.loc[eligible.NTC_detected_fraction.ge(.1)]
    def median_or_none(values):return float(values.median()) if len(values) else None
    result=dict(cells=len(cells),NTC_cells=int(cells.role.eq('non_targeting_control').sum()),
                designed_target_guides=len(mapping.loc[mapping.source_control_type.isin(['targeting','positive_control'])]),
                exact_id_targets=int(valid.target_gene_id.nunique()),missing_reference_targets=int(missing.target_gene_id.nunique()),
                guides_with_10cells=len(eligible),guides_NTC_detection_at_least_10pct=len(measured),
                median_RNA_remaining_ratio=median_or_none(eligible.RNA_remaining_ratio),
                median_ratio_NTC_detection_10pct=median_or_none(measured.RNA_remaining_ratio),
                RNA_not_DNA_editing=True,primary_batch_correction=False)
    atomic(out/'summary.json',result)
    atomic(out/'complete.json',dict(contract=contract,result=result,files={n:sha(out/n) for n in ['target_mapping.csv','guide_target_RNA.csv','summary.json','input_proof.json']}))
    return result


def run(config,resume=False):
    root=Path(config['output'])
    if root.parent!=Path(config['allowed_output_parent']):raise ValueError('Output parent differs')
    root.mkdir(exist_ok=resume)
    lock=root/'owner.lock'
    owner=dict(pid=os.getpid(),uid=os.getuid(),startticks=Path('/proc/self/stat').read_text().split(') ',1)[1].split()[19])
    with lock.open('x') as f:json.dump(owner,f)
    start=time.monotonic()
    def progress(stage,**kw):
        atomic(root/'status.json',dict(stage=stage,time=time.time(),elapsed=time.monotonic()-start,**owner,**kw));print(stage,kw,flush=True)
    def alarm(*_):raise TimeoutError('Time protection reached; completed data/checkpoints retained')
    signal.signal(signal.SIGALRM,alarm);signal.alarm(config['max_seconds'])
    try:
        results={}
        for d in config['datasets']:
            out=root/d['id'];complete=out/'complete.json'
            if complete.exists():
                r=json.loads(complete.read_text())
                if r['contract']!=hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest():raise ValueError('Dataset contract changed')
                for n,h in r['files'].items():
                    if sha(out/n)!=h:raise ValueError('Completed result changed')
                results[d['id']]=r['result'];continue
            progress('prepare_target_RNA',dataset=d['id'])
            results[d['id']]=run_dataset(d,out,lambda s,**kw:progress(s,dataset=d['id'],**kw))
        atomic(root/'summary.json',dict(datasets=results,code_commit=config['expected_commit'],no_DNA_editing_rate=True))
        progress('complete_pending_plot_review')
    except BaseException as exc:progress('failed_preserved',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    run(json.loads(a.config.read_text()),a.resume)
