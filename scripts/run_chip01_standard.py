"""Finite CPU first-pass QC/structure and exploratory perturbation statistics.

Read-only source counts; no arbitrary doublet cut, no invented replicates.
This is an adapted analysis, not an exact reproduction of either paper.
"""
from __future__ import annotations
import argparse, importlib.metadata as metadata, json, os, subprocess, time
from pathlib import Path
from datetime import datetime
import anndata as ad
import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
from analyze_expression_structure import blocks, sha256, atomic_json
from perturbation_statistics import differential_expression_block, bh_adjust, energy_v_statistic, energy_permutation_batched


def main_mask(obs):
    mask = obs.qc_pass_stage1.to_numpy(dtype=bool) & obs.guide_status.astype(str).eq('single_guide').to_numpy()
    if not np.array_equal(mask,obs.main_analysis_stage1.to_numpy(dtype=bool)):
        raise ValueError('Existing candidate flags disagree with explicit QC/single-guide policy')
    return mask


def normalize_counts(counts, totals, target_sum=10000):
    totals=np.asarray(totals,dtype=float)
    if np.any(totals<=0):raise ValueError('Nonpositive RNA counts')
    return (sparse.diags(target_sum/totals) @ counts.astype(np.float32)).tocsr().astype(np.float32)


def energy_representation(counts, totals, var, control, seed=0):
    """Mean >0.5 raw UMI, control-relative z, upper clip 10, 20 PCs.

    Adaptations: standard human ENSG/non-MT features, CP10k, no core-NTC
    selection. These choices and any zero-control-variance exclusions persist.
    """
    mean=np.asarray(counts.mean(axis=0)).ravel()
    genes=var.gene_id.astype(str).str.startswith('ENSG').to_numpy(copy=True)
    genes &= ~var.gene_name.astype(str).str.startswith('MT-').to_numpy()
    keep=np.flatnonzero(genes & (mean>.5))
    if len(keep)<21:raise ValueError('Insufficient expressed genes for energy representation')
    if len(totals)*len(keep)*4>2*1024**3:raise MemoryError('Selected-gene energy array exceeds 2GiB budget')
    normalized=normalize_counts(counts[:,keep],totals).toarray()
    center=normalized[control].mean(0);scale=normalized[control].std(0,ddof=1)
    valid=np.isfinite(scale)&(scale>1e-6)
    z=(normalized[:,valid]-center[valid])/scale[valid]
    np.minimum(z,10,out=z)
    pca=PCA(n_components=20,svd_solver='randomized',random_state=seed)
    scores=pca.fit_transform(z)
    return scores,keep[valid],pca.explained_variance_ratio_,int((~valid).sum())


def save_manifest(folder):
    manifest={p.relative_to(folder).as_posix():{'bytes':p.stat().st_size,'sha256':sha256(p)}
              for p in sorted(folder.rglob('*')) if p.is_file() and p.name not in ['status.json','output_manifest.json','owner.lock']}
    atomic_json(folder/'output_manifest.json',manifest)


def prepare(source,output,flagged_controls,diagnostics,state,seed):
    import scanpy as sc
    state('verify_input')
    expected='1dc39d7209faa23f8f099dd816f48fd59b3d1132dddea3cdd21accdafca4d8a0'
    if sha256(source)!=expected:raise ValueError('Input SHA256 mismatch')
    a=ad.read_h5ad(source,backed='r');obs=a.obs.copy();var=a.var.copy();shape=a.shape;a.file.close()
    if shape!=(25234,39663):raise ValueError('Unexpected source dimensions')
    selected=main_mask(obs);rows=np.flatnonzero(selected)
    parts=[]
    with h5py.File(source,'r') as f:
        if len(f['X/data'])!=122923364:raise ValueError('Source nnz mismatch')
        for picked,x in blocks(f,rows,len(var),128):
            if not np.array_equal(np.asarray(x.sum(1)).ravel(),obs.total_counts.to_numpy()[picked]):
                raise ValueError('Source row count mismatch')
            if (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():raise ValueError('Not nonnegative counts')
            parts.append(x)
    counts=sparse.vstack(parts,format='csr');del parts
    selected_obs=obs.iloc[rows].copy()
    flagged=set(pd.read_csv(flagged_controls).guide.astype(str))
    selected_obs['control_sensitivity_excluded']=selected_obs.is_control & selected_obs.guide_id.astype(str).isin(flagged)
    # Score annotation only; a missing score is not a confirmed singlet.
    score_table=pd.read_csv(diagnostics/'full_seed0/cell_doublet_diagnostics.csv.gz',index_col=0)
    score_col='doublet_score'
    selected_obs['doublet_diagnostic_score']=score_table[score_col].reindex(selected_obs.index)
    selected_obs['doublet_policy']='unresolved_score_no_additional_cut'
    selected_obs['analysis_policy']='rna_qc_single_guide_all_controls'
    selected_obs.to_csv(output/'cell_policy.csv.gz')
    coverage=selected_obs.groupby('perturbed_gene',observed=True).agg(cells=('guide_id','size'),guides=('guide_id','nunique'),control=('is_control','all'))
    coverage.to_csv(output/'coverage_per_target.csv')
    atomic_json(output/'qc_policy.json',{'source_sha256':expected,'original_cells':len(obs),'rna_qc_cells':int(obs.qc_pass_stage1.sum()),
      'main_cells':len(rows),'controls':int(selected_obs.is_control.sum()),'noncontrols':int((~selected_obs.is_control).sum()),
      'rna_thresholds':'reuse verified baseline qc_thresholds.json, not refitted',
      'single_guide_required':True,'extra_doublet_cut':None,'control_guides_removed_from_primary':[],
      'flagged_controls_for_sensitivity':sorted(flagged),'sensitivity_controls':int((selected_obs.is_control & ~selected_obs.control_sensitivity_excluded).sum()),
      'limitation':'Residual doublets possible; no support for existing automatic cutoff. Multi-guide exclusion is not complete doublet detection.',
      'source_modified':False,'biological_replicates':None,'dna_editing_efficiency':'not_assessed'})
    state('standard_structure',cells=len(rows))
    # Human expressed genes for representation; keep raw ALL features separately.
    human=var.gene_id.astype(str).str.startswith('ENSG').to_numpy(copy=True)
    human &= ~var.gene_name.astype(str).str.startswith('MT-').to_numpy()
    human &= np.asarray((counts>0).sum(0)).ravel()>=20
    struct=ad.AnnData(normalize_counts(counts[:,human],selected_obs.total_counts),obs=selected_obs.copy(),var=var.loc[human].copy())
    sc.pp.log1p(struct)
    sc.pp.highly_variable_genes(struct,n_top_genes=2000,flavor='seurat')
    sc.pp.pca(struct,n_comps=30,mask_var='highly_variable',svd_solver='arpack',random_state=seed)
    sc.pp.neighbors(struct,n_neighbors=15,n_pcs=30,random_state=seed)
    sc.tl.umap(struct,random_state=seed)
    sc.tl.leiden(struct,resolution=.5,flavor='igraph',n_iterations=2,directed=False,random_state=seed,key_added='expression_cluster')
    # Marker lists describe expression clusters, not cell-type annotation.
    sc.tl.rank_genes_groups(struct,'expression_cluster',method='wilcoxon',n_genes=30,use_raw=False)
    markers=sc.get.rank_genes_groups_df(struct,group=None)
    markers['gene_name']=markers.names.map(struct.var.gene_name)
    markers.to_csv(output/'expression_cluster_markers.csv.gz',index=False)
    coords=struct.obs[['perturbed_gene','guide_id','is_control','expression_cluster','total_counts','pct_counts_mt','doublet_diagnostic_score']].copy()
    coords[['UMAP1','UMAP2']]=struct.obsm['X_umap'];coords.to_csv(output/'cell_structure.csv.gz')
    pd.DataFrame(struct.obsm['X_pca'],index=selected_obs.index).to_csv(output/'expression_pca.csv.gz')
    struct.var.to_csv(output/'expression_features.csv')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(15,9),layout='constrained')
    for ax,key in zip(axes.ravel(),['expression_cluster','is_control','total_counts','pct_counts_mt','doublet_diagnostic_score','n_genes_by_counts']):
        vals=struct.obs[key]
        if key=='expression_cluster':vals=vals.cat.codes
        points=ax.scatter(*struct.obsm['X_umap'].T,c=vals,s=2,rasterized=True,cmap='viridis')
        ax.set_title(key);ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');fig.colorbar(points,ax=ax,shrink=.7)
    fig.savefig(output/'expression_umap.png',dpi=160);plt.close(fig)
    # Compact graph/embedding object: avoids another full duplicated expression matrix.
    graph=ad.AnnData(sparse.csr_matrix((len(rows),0)),obs=struct.obs.copy())
    for key in ['X_pca','X_umap']:graph.obsm[key]=struct.obsm[key]
    graph.obsp['connectivities']=struct.obsp['connectivities'];graph.obsp['distances']=struct.obsp['distances']
    graph.uns['counts_source']=str(source);graph.uns['counts_source_sha256']=expected
    graph.write_h5ad(output/'cell_structure.h5ad',compression='gzip')
    del struct
    state('formal_energy_representation',cells=len(rows))
    control=selected_obs.is_control.to_numpy(dtype=bool)
    pcs,features,variance,excluded=energy_representation(counts,selected_obs.total_counts,var,control,seed)
    np.savez_compressed(output/'energy_representation.npz',coordinates=pcs,feature_indices=features,variance_ratio=variance,source_rows=rows)
    var.iloc[features].to_csv(output/'energy_features.csv')
    atomic_json(output/'representation_parameters.json',{'expression':'log1p CP10k, ENSG non-MT expressed>=20 cells; Seurat 2000 HVG; PCA30; neighbors15; UMAP; Leiden0.5',
      'energy':'ENSG non-MT mean raw UMI>0.5; CP10k; z relative to all NTC; upper clip10; PCA20 randomized',
      'energy_genes':len(features),'zero_control_variance_genes_excluded':excluded,
      'adaptations_vs_Replogle2022':['human ENSG/non-MT only','CP10k depth target','all NTC retained, no core-NTC selection','bounded permutation subsampling in effects stage'],
      'exact_paper_replication':False,'seed':seed})


def effects(source,output,prepared,state,seed,permutations):
    cp=np.load(prepared/'energy_representation.npz');rows=cp['source_rows'];pcs=cp['coordinates']
    obs=pd.read_csv(prepared/'cell_policy.csv.gz',index_col=0)
    a=ad.read_h5ad(source,backed='r');var=a.var.copy();a.file.close()
    with h5py.File(source,'r') as f:counts=sparse.vstack([x for _,x in blocks(f,rows,len(var))],format='csr')
    if not np.array_equal(np.asarray(counts.sum(1)).ravel(),obs.total_counts.to_numpy()):raise ValueError('Count totals changed')
    human=var.gene_id.astype(str).str.startswith('ENSG').to_numpy()
    # Common tested universe, not a different filter chosen for each target.
    features=np.flatnonzero(human & (np.asarray((counts>0).sum(0)).ravel()>=20))
    normalized=normalize_counts(counts[:,features],obs.total_counts)
    del counts
    control=obs.is_control.to_numpy(dtype=bool)
    alternative=control & ~obs.control_sensitivity_excluded.to_numpy(dtype=bool)
    labels=obs.perturbed_gene.astype(str).to_numpy()
    sizes=pd.Series(labels[~control]).value_counts();targets=sorted(sizes[sizes>=10].index)
    var.iloc[features].to_csv(output/'tested_gene_universe.csv')
    ctrl=normalized[control].tocsc();alt=normalized[alternative].tocsc()
    rng=np.random.default_rng(seed)
    reference=np.sort(rng.choice(np.flatnonzero(control),min(256,control.sum()),replace=False))
    alt_reference=np.sort(rng.choice(np.flatnonzero(alternative),min(256,alternative.sum()),replace=False))
    np.savez_compressed(output/'energy_reference_rows.npz',primary=reference,sensitivity=alt_reference)
    within_control=float(cdist(pcs[control],pcs[control]).mean())
    within_alternative=float(cdist(pcs[alternative],pcs[alternative]).mean())
    energies=[];summary=[];profiles=[];started=time.monotonic()
    for i,target in enumerate(targets):
        state('formal_effects',completed_targets=i,total_targets=len(targets),target=target,elapsed_seconds=time.monotonic()-started)
        idx=np.flatnonzero((labels==target)&~control);x=normalized[idx].tocsc();tables=[]
        for start in range(0,len(features),128):
            stop=min(start+128,len(features));dense=x[:,start:stop].toarray()
            primary=differential_expression_block(dense,ctrl[:,start:stop].toarray())
            sensitivity=differential_expression_block(dense,alt[:,start:stop].toarray())
            table=pd.DataFrame(primary,index=var.index[features[start:stop]])
            for key in ['pvalue','log2_fold_change']:table['ntc_sensitivity_'+key]=sensitivity[key]
            tables.append(table)
        table=pd.concat(tables);table['bh_q']=bh_adjust(table.pvalue.to_numpy())
        table['ntc_sensitivity_bh_q']=bh_adjust(table.ntc_sensitivity_pvalue.to_numpy())
        table['gene_name']=var.gene_name.reindex(table.index);table.index.name='gene_id'
        significant=(table.bh_q<.05)&(table.log2_fold_change.abs()>=.25)
        table['exploratory_de']=significant
        name=f'target_{i:03d}'
        table.to_csv(output/(name+'_de.csv.gz'))
        profiles.append(table.log2_fold_change.to_numpy())
        summary.append(dict(file=name+'_de.csv.gz',target=target,cells=len(idx),control_cells=int(control.sum()),
                            tested_genes=len(features),exploratory_de_genes=int(significant.sum()),
                            sensitivity_de_genes=int(((table.ntc_sensitivity_bh_q<.05)&(table.ntc_sensitivity_log2_fold_change.abs()>=.25)).sum()),
                            low_coverage=len(idx)<20))
        if len(idx)>=20:
            chosen=np.sort(rng.choice(idx,min(256,len(idx)),replace=False))
            # Full-data descriptive effect plus bounded permutation p-value.
            within_target=float(cdist(pcs[idx],pcs[idx]).mean())
            descriptive=max(0.0,float(2*cdist(pcs[idx],pcs[control]).mean()-within_target-within_control))
            sampled,p=energy_permutation_batched(pcs[chosen],pcs[reference],permutations=permutations,seed=seed+i)
            sensitivity=max(0.0,float(2*cdist(pcs[idx],pcs[alternative]).mean()-within_target-within_alternative))
            sampled_alt,p_alt=energy_permutation_batched(pcs[chosen],pcs[alt_reference],permutations=permutations,seed=seed+i)
            energies.append(dict(target=target,cells=len(idx),sampled_target_cells=len(chosen),sampled_control_cells=len(reference),
                energy_full=descriptive,energy_sampled=sampled,pvalue=p,ntc_sensitivity_energy_full=sensitivity,
                ntc_sensitivity_energy_sampled=sampled_alt,ntc_sensitivity_pvalue=p_alt,permutations=permutations))
        pd.DataFrame(summary).to_csv(output/'de_summary.csv',index=False)
        pd.DataFrame(energies).to_csv(output/'energy_progress.csv',index=False)
        atomic_json(output/'progress.json',dict(completed_targets=i+1,total_targets=len(targets),elapsed_seconds=time.monotonic()-started,
             seconds_per_target=(time.monotonic()-started)/(i+1)))
    energy=pd.DataFrame(energies)
    if len(energy):
        energy['bh_q']=bh_adjust(energy.pvalue.to_numpy());energy['ntc_sensitivity_bh_q']=bh_adjust(energy.ntc_sensitivity_pvalue.to_numpy())
    energy.to_csv(output/'energy_results.csv',index=False)
    np.savez_compressed(output/'perturbation_expression_profiles.npz',values=np.array(profiles),targets=np.array(targets),gene_indices=features)
    atomic_json(output/'statistics_parameters.json',{'de':'two-sided asymptotic tie-corrected Mann-Whitney; per-target BH over entire common gene universe',
       'de_min_cells':10,'gene_universe':'human ENSG detected in>=20 primary cells, all such genes tested',
       'de_effect':'CP10k mean log2 fold change with pseudocount1; display q<.05 and abs(log2FC)>=.25',
       'energy_min_cells':20,'energy_descriptive':'all target cells vs all primary/alternative controls',
       'energy_test':'fixed shared seeded controls<=256; each target<=256; no replacement; batched energy V permutation',
       'permutations':permutations,'seed':seed,'biological_replication':False,'inference_scope':'exploratory within-chip cell-level association, not replicated causal inference',
       'sensitivity':'exclude 8 previously flagged control guides only; same target cells and representation',
       'sensitivity_limitation':'NTC-exclusion sensitivity does not establish robustness to unresolved doublets'})


def run(args):
    if args.output.exists():raise FileExistsError('Inspect existing stage; never overwrite/restart blindly')
    if psutil.virtual_memory().available<16*1024**3:raise MemoryError('Need 16GiB available on execution host')
    args.output.mkdir(parents=True)
    lock=args.output/'owner.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    started=time.monotonic();commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    def state(stage,**kw):atomic_json(args.output/'status.json',dict(stage=stage,pid=os.getpid(),time=datetime.now().astimezone().isoformat(),code_commit=commit,**kw))
    try:
        with threadpool_limits(2):
            if args.stage=='prepare':prepare(args.source,args.output,args.flagged_controls,args.diagnostics,state,args.seed)
            else:effects(args.source,args.output,args.prepared,state,args.seed,args.permutations)
        save_manifest(args.output)
        state('complete',elapsed_seconds=time.monotonic()-started,source_modified=False,
              versions={k:metadata.version(k) for k in ['scanpy','anndata','numpy','scipy','scikit-learn']},
              next_step='effects' if args.stage=='prepare' else 'functional_analysis_and_report')
    except Exception as exc:state('failed',error=repr(exc));raise
    finally:lock.unlink()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['prepare','effects'],required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--flagged-controls',type=Path);p.add_argument('--diagnostics',type=Path);p.add_argument('--prepared',type=Path)
    p.add_argument('--seed',type=int,default=0);p.add_argument('--permutations',type=int,default=9999)
    run(p.parse_args())
