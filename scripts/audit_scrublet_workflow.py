"""Transparent native Scanpy/Scrublet preprocessing audit, never cell filtering.

Starts from RNA-QC-passing raw counts, lets the native HVG procedure select
features once, and records that selection. Scores and automatic calls remain
diagnostic outputs, not a validated main-analysis mask.
"""
from __future__ import annotations
import argparse,json,os,subprocess,time
from pathlib import Path
import importlib.metadata as metadata
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from sklearn import config_context
from threadpoolctl import threadpool_limits
from analyze_expression_structure import atomic_json,blocks,sha256


def prepare_native_inputs(counts,obs,var,seed=0,sim_ratio=2.0):
    """Match Scanpy 1.12's no-adata_sim preparation using public functions.

    Keeping this explicit exposes filtered/selected dimensions. Passing the
    resulting adata_sim to scrublet avoids selecting HVGs a second time.
    """
    import scanpy as sc
    if not sparse.issparse(counts):
        raise ValueError('Raw input must be sparse')
    if not np.isfinite(counts.data).all() or (counts.data<0).any() or not np.equal(counts.data,np.floor(counts.data)).all():
        raise ValueError('Expected finite nonnegative integer RNA counts')
    data=ad.AnnData(counts.copy(),obs=obs.copy(),var=var.copy())
    before=list(data.shape)
    sc.pp.filter_genes(data,min_cells=3)
    sc.pp.filter_cells(data,min_genes=3)
    filtered=list(data.shape)
    if len(data)!=len(obs):
        raise ValueError('Native preprocessing excluded input cells; inspect instead of silently changing the QC set')
    data.layers['raw']=data.X.copy()
    sc.pp.normalize_total(data)
    data.layers['log1p']=data.X.copy()
    sc.pp.log1p(data,layer='log1p')
    sc.pp.highly_variable_genes(data,layer='log1p')
    del data.layers['log1p']
    data=data[:,data.var.highly_variable].copy()
    simulated=sc.pp.scrublet_simulate_doublets(data,layer='raw',sim_doublet_ratio=sim_ratio,random_seed=seed)
    del data.layers['raw']
    sc.pp.normalize_total(data,target_sum=1e6)
    sc.pp.normalize_total(simulated,target_sum=1e6)
    info={'input_shape':before,'after_native_minimum_filters':filtered,
          'selected_hvg_count':data.n_vars,'simulated_cells':simulated.n_obs,
          'feature_selection_uses_guide_labels':False,
          'recipe':'Scanpy native: min_cells3/min_genes3; save raw; normalize; log1p/HVG; raw-count simulation; normalize both to1e6; advanced scrublet with supplied simulations',
          'preselected_diagnostic_features':False,'hvg_selection_passes':1}
    return data,simulated,info


def run(source,output,expected_sha256,seed=0,limit=0):
    if output.exists():
        raise FileExistsError('Preserve existing attempt; inspect before choosing a new output')
    output.mkdir(parents=True)
    fd=os.open(output/'owner.lock',os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    started=time.monotonic()
    def state(stage,**kw):
        atomic_json(output/'status.json',dict(stage=stage,pid=os.getpid(),elapsed_seconds=time.monotonic()-started,**kw))
    try:
        state('verify_raw_source')
        if sha256(source)!=expected_sha256:raise ValueError('Source hash mismatch')
        if psutil.virtual_memory().available< (4 if limit else 16)*1024**3:
            raise MemoryError('Need 4GiB smoke /16GiB full available RAM before native all-RNA audit')
        a=ad.read_h5ad(source,backed='r');obs=a.obs.copy();var=a.var.copy();n_vars=a.n_vars;a.file.close()
        rows=np.flatnonzero(obs.qc_pass_stage1.to_numpy(dtype=bool))
        if limit:rows=np.sort(np.random.default_rng(0).choice(rows,min(limit,len(rows)),replace=False))
        state('read_raw_counts',cells=len(rows),rna_features=n_vars)
        with h5py.File(source,'r') as f:
            parts=[x for _,x in blocks(f,rows,n_vars)]
        counts=sparse.vstack(parts,format='csr');del parts
        assert np.array_equal(np.asarray(counts.sum(1)).ravel(),obs.total_counts.to_numpy()[rows])
        state('native_preprocessing',cells=len(rows),rna_features=n_vars)
        with config_context(working_memory=64),threadpool_limits(limits=2):
            data,sim,info=prepare_native_inputs(counts,obs.iloc[rows],var,seed)
            del counts
            if data.n_vars<=30:raise ValueError('Too few native HVGs for30PC diagnostic; no silent parameter change')
            atomic_json(output/'preprocessing.json',info)
            data.var.to_csv(output/'selected_features.csv')
            state('score_scrublet',**info)
            import scanpy as sc
            sc.pp.scrublet(data,adata_sim=sim,expected_doublet_rate=.05,sim_doublet_ratio=2.,mean_center=True,n_prin_comps=30,use_approx_neighbors=False,random_state=seed)
        result=data.uns['scrublet'];scores=data.obs.doublet_score.to_numpy();synthetic=np.asarray(result['doublet_scores_sim'])
        assert np.isfinite(scores).all() and np.isfinite(synthetic).all()
        data.obs.to_csv(output/'cell_doublet_diagnostics.csv.gz')
        np.save(output/'simulated_doublet_scores.npy',synthetic)
        threshold=result.get('threshold')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4))
        ax.hist(scores,bins=60,density=True,alpha=.6,label='Observed records; identities unknown')
        ax.hist(synthetic,bins=60,density=True,alpha=.6,label='Synthetic doublets')
        if threshold is not None:ax.axvline(threshold,color='black',ls='--',label='Automatic candidate cutoff; not validated')
        ax.set(xlabel='Diagnostic doublet score',ylabel='Density',title='Native all-RNA input audit; no cell exclusion')
        ax.legend(fontsize=8);fig.tight_layout();fig.savefig(output/'score_distribution.png',dpi=160);plt.close(fig)
        state('complete_smoke_only' if limit else 'complete_diagnostics_only',**info,
              cells=len(rows),seed=seed,expected_rate_assumption=.05,experimental_doublet_rate='unverified',
              candidate_threshold=float(threshold) if threshold is not None else None,
              candidate_calls=int((scores>threshold).sum()) if threshold is not None else None,
              simulated_fraction_above_threshold=float((synthetic>threshold).mean()) if threshold is not None else None,
              input_sha256=expected_sha256,code_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
              versions={p:metadata.version(p) for p in ['scanpy','scikit-image','numpy','scipy','anndata']},
              filtering_applied=False,source_counts_modified=False,final_qc=False,
              output_sha256={p.name:sha256(p) for p in output.iterdir() if p.name!='owner.lock' and p.name!='status.json'})
    except Exception as e:state('failed',error=repr(e));raise
    finally:(output/'owner.lock').unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--expected-sha256',required=True);p.add_argument('--seed',type=int,default=0);p.add_argument('--limit',type=int,default=0)
    a=p.parse_args();run(a.source,a.output,a.expected_sha256,a.seed,a.limit)
