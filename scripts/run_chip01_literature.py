"""Staged CPU literature-aligned redo, immutable source and distinct outputs."""
from pathlib import Path
from datetime import datetime
import argparse,csv,hashlib,json,os,subprocess,time,importlib.metadata as md
import numpy as np
import pandas as pd
import anndata as ad
import h5py,psutil
from scipy import sparse
from threadpoolctl import threadpool_limits
from analyze_expression_structure import blocks,sha256,atomic_json
from perturbation_statistics import bh_adjust
from replogle2022_statistics import ADCalibration,ad_statistics


def read_input(source,expected):
    if sha256(source)!=expected:raise ValueError('Source SHA256 mismatch')
    a=ad.read_h5ad(source,backed='r');obs=a.obs.copy();var=a.var.copy();shape=a.shape;a.file.close()
    selected=obs.qc_pass_stage1.astype(bool)&obs.guide_status.astype(str).eq('single_guide')
    if not np.array_equal(selected,obs.main_analysis_stage1.astype(bool)):raise ValueError('Existing candidate policy mismatch')
    rows=np.flatnonzero(selected.to_numpy())
    with h5py.File(source,'r') as f:
        if f['X'].attrs['encoding-type']!='csr_matrix':raise ValueError('Sparse CSR source required')
        counts=sparse.vstack([x for _,x in blocks(f,rows,len(var))],format='csr')
        nnz=len(f['X/data'])
    meta=obs.iloc[rows].copy();meta['source_row']=rows
    if not np.array_equal(np.asarray(counts.sum(1)).ravel(),meta.total_counts.to_numpy()):raise ValueError('Input count invariants differ')
    if (counts.data<0).any() or np.any(counts.data!=np.floor(counts.data)):raise ValueError('Expected original integer counts')
    return counts,meta,var,dict(source_shape=list(shape),source_nnz=nnz,selected_cells=len(rows),sha256=expected,source_modified=False)


def validate_calibration(calibration,validation):
    table=pd.read_csv(validation)
    observed=calibration.pvalues(table.statistic)
    error=np.max(np.abs(observed-table.pvalue.to_numpy())/np.maximum(table.pvalue.to_numpy(),1e-300))
    if error>.001:raise ValueError(f'R AD off-grid relative error {error} exceeds 0.1% gate')
    return dict(reference_points=len(table),maximum_relative_error=float(error),calibration_validated=True)


def core_controls(args,config,state):
    state('verify_calibration')
    calibration=ADCalibration(args.calibration)
    validation=validate_calibration(calibration,args.calibration_validation)
    atomic_json(args.output/'calibration_validation.json',validation)
    state('verify_and_read_sparse_source')
    counts,obs,var,identity=read_input(args.source,args.input_sha256)
    control=obs.is_control.to_numpy(dtype=bool)
    guides=obs.guide_id.astype(str).to_numpy()
    ctrl=counts[control]
    means=np.asarray(ctrl.mean(0)).ravel()
    features=np.flatnonzero(means>config['core_control_gene_mean_raw_umi_gt'])
    if not len(features):raise ValueError('No genes satisfy paper core-control expression threshold')
    sizes=pd.Series(guides[control]).value_counts()
    candidates=sorted(sizes[sizes>=config['core_control_min_cells_for_pairwise_ad']].index)
    if len(candidates)<2:raise ValueError('Too few mathematically testable core-control candidates')
    sorted_groups={g:np.sort(counts[control&(guides==g)][:,features].toarray().astype(float),axis=0) for g in candidates}
    total=len(candidates)*(len(candidates)-1)//2;completed=0;started=time.monotonic()
    sums={g:0 for g in candidates};pair_counts={g:0 for g in candidates}
    with (args.output/'control_pairwise_ad.csv').open('w',encoding='utf8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['guide_a','guide_b','cells_a','cells_b','tested_genes','ad_bh_de_genes']);writer.writeheader()
        for i,left in enumerate(candidates):
            for right in candidates[i+1:]:
                statistics=ad_statistics(sorted_groups[left],sorted_groups[right],already_sorted=True)
                p=calibration.pvalues(statistics);q=bh_adjust(p);de=int((q<config['gene_bh_q']).sum())
                writer.writerow(dict(guide_a=left,guide_b=right,cells_a=len(sorted_groups[left]),cells_b=len(sorted_groups[right]),tested_genes=len(features),ad_bh_de_genes=de))
                sums[left]+=de;sums[right]+=de;pair_counts[left]+=1;pair_counts[right]+=1;completed+=1
                if completed%16==0 or completed==total:
                    stream.flush()
                    state('core_control_pairwise_ad',completed_pairs=completed,total_pairs=total,
                          tested_genes=len(features),elapsed_seconds=time.monotonic()-started,
                          seconds_per_pair=(time.monotonic()-started)/completed)
    records=[]
    for guide,size in sizes.items():
        testable=guide in pair_counts
        avg=sums[guide]/pair_counts[guide] if testable else None
        keep=testable and avg<config['core_control_mean_ad_de_genes_lt']
        records.append(dict(guide_id=guide,cells=int(size),pairwise_ad_testable=testable,mean_ad_de_genes=avg,core_control=keep,
                            reason='below_paper_core_threshold' if keep else ('above_paper_core_threshold' if testable else 'insufficient_cells_for_pairwise_AD')))
    summary=pd.DataFrame(records).sort_values('guide_id')
    summary.to_csv(args.output/'core_control_guides.csv',index=False)
    core_guides=set(summary.loc[summary.core_control,'guide_id'])
    obs['core_control']=control&np.isin(guides,list(core_guides))
    if obs.core_control.sum()<25:raise ValueError('Fewer than 25 core control cells; inspect results, do not change threshold automatically')
    obs['literature_analysis_cell']=~control|obs.core_control.to_numpy()
    obs.to_csv(args.output/'cell_policy.csv.gz')
    var.iloc[features].to_csv(args.output/'core_control_tested_genes.csv',index=False)
    identity.update(core_control_guides=len(core_guides),original_control_guides=len(sizes),core_control_cells=int(obs.core_control.sum()),
                    candidate_guides=len(candidates),testable_pairs=total,core_threshold=config['core_control_mean_ad_de_genes_lt'],
                    normalization_target_umi=float(np.median(obs.loc[obs.core_control,'total_counts'])),
                    calibration_sha256=sha256(args.calibration),calibration_validation=validation)
    atomic_json(args.output/'core_summary.json',identity)


def run(args):
    if args.output.exists():raise FileExistsError('Distinct run output required; inspect before recovery')
    if psutil.virtual_memory().available<12*1024**3:raise MemoryError('Execution host needs 12GiB free; do not run full counts locally')
    config=json.loads(args.config.read_text())
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    start=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**kw):atomic_json(args.output/'status.json',dict(stage=stage,pid=os.getpid(),time=datetime.now().astimezone().isoformat(),code_commit=commit,**kw))
    try:
        atomic_json(args.output/'frozen_parameters.json',config)
        with threadpool_limits(config['threads']):
            if args.stage=='core':core_controls(args,config,state)
            else:raise ValueError('Unsupported stage')
        atomic_json(args.output/'output_manifest.json',{p.name:dict(bytes=p.stat().st_size,sha256=sha256(p)) for p in args.output.iterdir() if p.is_file() and p.name not in ['status.json','owner.lock','output_manifest.json']})
        state('complete',elapsed_seconds=time.monotonic()-start,stage_completed=args.stage,
              versions={k:md.version(k) for k in ['numpy','scipy','numba','anndata','scikit-learn']},source_modified=False)
    except Exception as error:state('failed',error=repr(error));raise
    finally:lock.unlink()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['core'],required=True)
    for name in ['source','output','config','calibration','calibration-validation']:parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--input-sha256',required=True)
    run(parser.parse_args())
