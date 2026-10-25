"""Bounded comparison of raw local signatures and NTC-relative contrasts.

This audits a representation mismatch. It is not a full Mixscape reproduction,
KO classification, a new discovery screen or a negative-control calibration.
"""
from pathlib import Path
import argparse,gc,json,os,signal,subprocess,time
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

from compare_established_perturbseq import load_counts,ntc_z_parameters,nearest_controls,group_sets,evaluate


def run(config_path):
 config=json.loads(Path(config_path).read_text());out=Path(config['output']);out.mkdir(parents=True,exist_ok=False)
 lock=out/'owner.lock';lock.write_text(str(os.getpid()));start=time.monotonic()
 def status(stage,**kwargs):
  value=dict(stage=stage,pid=os.getpid(),elapsed_seconds=time.monotonic()-start,**kwargs)
  p=out/'status.tmp';p.write_text(json.dumps(value));p.replace(out/'status.json');print(json.dumps(value),flush=True)
 def deadline(signum,frame):raise TimeoutError('Local-signature audit deadline')
 signal.signal(signal.SIGALRM,deadline);signal.alarm(config.get('max_seconds',360))
 try:
  with threadpool_limits(limits=2):
   status('load_count_identical_input')
   x,cells,genes,proof=load_counts(config,status);(out/'input_proof.json').write_text(json.dumps(proof,indent=2))
   ntc=cells.role.eq('non_targeting_control').to_numpy();totals=np.asarray(x.sum(1)).ravel();scale=float(np.median(totals[ntc]))
   norm=x.astype(np.float32,copy=False);del x
   norm.data*=np.repeat((scale/totals).astype(np.float32),np.diff(norm.indptr))
   _,sd,valid=ntc_z_parameters(norm,ntc)
   norm.data=np.log1p(norm.data*(10000/scale)).astype(np.float32)
   a=ad.AnnData(norm,var=pd.DataFrame(index=genes.gene_id))
   sc.pp.highly_variable_genes(a,flavor='seurat',n_top_genes=2000,inplace=True)
   hvg=np.flatnonzero(a.var.highly_variable.to_numpy()&valid);matrix=norm[:,hvg].toarray();del a,norm;gc.collect()
   if 'expected_hvg' in config:assert len(hvg)==config['expected_hvg']
   status('local_signature',hvg=len(hvg))
   embedding=PCA(n_components=40,svd_solver='randomized',random_state=20261026).fit_transform(matrix)
   neighbors=nearest_controls(embedding,ntc,cells.guide_id.to_numpy(),k=20)
   residual=matrix.copy()
   for first in range(0,len(cells),256):residual[first:first+256]-=matrix[neighbors[first:first+256]].mean(1)
   control_mean=residual[ntc].mean(0)
   pd.DataFrame(dict(gene_id=genes.gene_id.iloc[hvg].to_numpy(),ntc_mean_local_signature=control_mean)).to_csv(out/'ntc_mean_local_signature.csv.gz',index=False)
   groups=group_sets(cells,20261026);results=[];saved={}
   for name,mode in [('raw_local_signature','signature_mean'),('ntc_relative_local_contrast','difference')]:
    status('evaluate',quantity=name)
    summary,tables,effects=evaluate(name,residual,cells,genes.iloc[hvg],groups,mode,None,hvg,return_guide_profiles=True)
    saved[name]=dict(table=tables['guide_groups'],effects=effects)
    np.savez_compressed(out/(name+'_guide_profiles.npz'),effects=effects,groups=tables['guide_groups'].group.to_numpy(dtype=str),gene_ids=genes.gene_id.iloc[hvg].to_numpy(dtype=str),ntc_mean_signature=control_mean)
    summary['group_profile_rms_median']=float(np.median(np.sqrt(np.mean(effects**2,axis=1))))
    results.append(summary)
    for key,table in tables.items():table.to_csv(out/(name+'_'+key+'.csv.gz'),index=False)
   # The raw quantity must reproduce its previous implementation, not merely
   # look similar after changes to centering or background-role selection.
   if 'prior_same_target_pairs' in config:
    prior=pd.read_csv(config['prior_same_target_pairs']);current=pd.read_csv(out/'raw_local_signature_same_target_pairs.csv.gz')
    joined=current.merge(prior[['a','b','r']],on=['a','b'],suffixes=('_now','_prior'),validate='one_to_one')
    assert len(joined)==len(prior)==len(current)
    maximum=float(np.max(np.abs(joined.r_now-joined.r_prior)))
    assert maximum<1e-6,'Raw-signature reproduction differs'
   else:maximum=None
   final=dict(code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[1],text=True).strip(),
    config=config,cells=len(cells),ntc_cells=int(ntc.sum()),genes=len(genes),hvg=len(hvg),
    normalization_target_umi=scale,ntc_signature_mean_rms=float(np.sqrt(np.mean(control_mean**2))),
    raw_signature_prior_max_r_difference=maximum,representations=results,
    max_rss_kib=__import__('resource').getrusage(__import__('resource').RUSAGE_SELF).ru_maxrss,
    elapsed_seconds=time.monotonic()-start,interpretation='Quantity-comparability audit only; no method-power or biological success claim')
   (out/'summary.json').write_text(json.dumps(final,indent=2));status('complete')
 except BaseException as error:status('failed',error=repr(error));raise
 finally:signal.alarm(0);lock.unlink()


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--config',required=True);args=p.parse_args();run(args.config)
